## Planned improvements (not bugs, just known follow-up work)

- **`Problem` type decoupling mesh/point generation from basis — done, as a thin
  wrapper.** `MechanicalProblem` bundles `Mesh`+`Point`+`Time` — the expensive,
  IC-defining part of a simulation — as one object, reusing the existing
  `MaterialPointPhase` hierarchy rather than a second "phase" concept. `Basis`/`Solver`
  remain separate, un-bundled types, since `basis.which` and `basis.trsfr` vary
  independently. `types/problem/` holds every file that feeds `MechanicalProblem`
  construction — `constitutive.jl`, `geometry.jl`, `eulerian.jl` (Mesh), `tensor.jl`
  (the `AbstractStrain`/`AbstractStress` tensor types `Point` stores), `lagrangian.jl`
  (Point), `time.jl` (Time), `problem.jl` (the bundle itself) — loaded in that explicit
  order by `boot.jl`'s hand-ordered include list (see the main `CLAUDE.md`'s "Project
  shape" section), not by alphabetical file order.

  `setup_problem` (`src/home/init/setup_problem.jl`) is the entry point — calls
  `setup_mesh`/`setup_material_constants`/`setup_mpts`/`setup_time` internally and
  returns a `Problem`; `slump_problem` (renamed from `ic_slump`) uses it.
  `export_problem` writes a single `ic["problem"]` entry instead of separate
  `ic["mesh"]`/`ic["mpts"]`/`ic["time"]` keys — see `persistence.md`.
  `Basis` stays its own top-level `ic["basis"]` key, deliberately independent of
  `Problem`. Kernel/workflow signatures (`workflow!(mpts,mesh,basis,time,solver)`) are
  untouched — `Problem` is, for now, only an outer-API/persistence convenience, not
  threaded through kernels. `setup_basis` takes `(problem::MechanicalProblem, solver)`.
- **Kernel granularity**: `elast!` is the model to follow — decomposed into small
  functions each doing one specific task, rather than one monolithic kernel body.
  `update.jl`'s `elastoplast`/`elasto` dispatch functions (branch on
  `material.elastic`/`basis.how` inline) would benefit from the same split.
- **Typed constitutive-model abstraction — done.** `mpts.s.cmp::Vector{CM} where
  CM<:AbstractConstitutiveModel` (`constitutive.jl`) bundles the *static* elastic+
  plastic material constants (`Gc`, `Kc`, `Del`, `Hp`, `c₀`, `cᵣ`, `ϕ₀`) into one typed
  object per particle, replacing an old global `cmpr::NamedTuple` threaded through ~15
  signatures. Evolving plastic state (`Δλ`, `ϵpII`, `ϵpV`) deliberately stays *outside*
  `cmp`, as flat mutable vectors on `PointSolidPhase`, since it's written every
  timestep and bundling it into an immutable per-particle struct would force a
  full-struct reconstruction on every write.
  - `PerfectlyElastic` is defined but has **no construction path** — dead scaffolding.
    `DruckerPrager`/`VonMises` both **do** have construction paths, branched on
    `material.plastic` — see "DP/J2 retmap kernel unification" below.
  - `PointSolidPhase.rheo::R`/`AbstractRheology` — a parallel, never-read earlier
    attempt at the same per-particle constitutive-data problem — was retired entirely,
    since `cmp` already supersedes it.
  - `dynamic_relaxation`'s `cmp`-rerouted reads are confirmed working for the plain
    (non-u-P) path via `test_column.jl`'s convergence sweep, including its
    correctness assertions. The u-P variant remains unverified (and throws immediately
    if exercised regardless — see `.claude/bug/known/dynamic-relaxation-uP-fint-p2n-missing.md`).
- **DP/J2 retmap kernel unification — done.** `retmap/DP.jl`/`retmap/J2.jl`'s four
  separately-named kernels (`finite_DP`/`infinitesimal_DP`/`finite_J2`/
  `infinitesimal_J2`, picked between via a `plast.constitutive`/`strain.deform` string
  branch) are now one kernel name, `retmap`, with four methods dispatched on `Point`'s
  `CM` (`DruckerPrager`/`VonMises`) and `ST` (`LogarithmicStrain`/`InfinitesimalStrain`)
  — mirroring `elast!`'s `ST`-only dispatch, extended to a second axis.
  `init_update` registers `retmap(CPU())` unconditionally, no string branch.

  `VonMises{T2,D,NSTR,L} <: AbstractConstitutiveModel{T2,D}` — same fields as
  `DruckerPrager` minus `ϕ₀` (J2 yield has no pressure dependence). `setup_cmp`
  branches on a `constitutive::String` keyword to build `Vector{DruckerPrager}`
  (`"DP"`) or `Vector{VonMises}` (`"VM"`); an unrecognized string throws immediately at
  setup time. Config string unified to `"VM"` everywhere. Both `DP.jl`/`J2.jl` share a
  per-model return-mapping core (`_druckerprager_return_map`/`_vonmises_return_map`)
  called by both the finite- and infinitesimal-strain `drucker_prager`/`von_mises`
  methods, rather than the infinitesimal kernels hand-duplicating the algebra
  separately — verified the extracted CPA loop matches the original inline version
  bit-for-bit.

  **Gotcha**: `slump_problem` only plots `"phi0"` (friction angle) in its initial-
  condition figure when `material.plastic=="DP"` — `VonMises` genuinely has no `ϕ₀`
  field, so plotting it unconditionally used to crash under `"VM"`. Don't fake a value;
  the field just isn't applicable for a pressure-independent yield surface.
- **Typed strain/stress tensor storage — done.** `mpts.s.ϵᵢⱼ`/`.ϵn`/`.σᵢⱼ`/`.σn`/`.τᵢⱼ`
  hold typed tensor objects from `src/boot/needs/types/problem/tensor.jl`, each
  storing a volumetric+deviatoric additive split:
  - `LogarithmicStrain{S,T,L}` / `InfinitesimalStrain{S,T,L}` (`vol::T`,
    `dev::SMatrix{S,S,T,L}`) — `ϵᵢⱼ`/`ϵn`, picked by `solver.material.elastic`.
  - `CauchyStress{S,T,L}` / `KirchhoffStress{S,T,L}` (`p::T` positive in compression,
    `dev::SMatrix{S,S,T,L}`) — `σᵢⱼ`/`σn` are always Cauchy, `τᵢⱼ` always Kirchhoff.

  The abstract supertypes (`AbstractTensor`/`AbstractStrain`/`AbstractStress`) live in
  `types/problem/tensor.jl` alongside their concrete subtypes — `tensor.jl` is
  deliberately included before `types/problem/lagrangian.jl` in `boot.jl`'s explicit
  list (see the main `CLAUDE.md`'s "Project shape" section), since `Point`/
  `PointSolidPhase` there constrain `ST` against `AbstractStrain` and spell their stress
  fields as `CauchyStress{D,T2,L}`/`KirchhoffStress{D,T2,L}` (no separate `SC`/`SK`
  parameters — they never varied).

  Two things worth knowing before touching this:
  - **The stored split is representational, not canonical.** The only invariant
    guaranteed is that `get_voigt`/`get_tensor` reconstruct the tensor. `dev` is *not*
    required to be trace-free, because two different pressure conventions coexist
    (`_trial_elastic_stress` uses `tr(ϵ)/3` even in 2D plane strain; the Voigt
    constructors use `-tr(σ)/S`, matching `σTr`/`yield_J2`'s divide-by-2-in-2D).
    `get_J2`/`get_τII` re-derive a trace-free deviator from `get_voigt` rather than
    trusting the stored `dev`. **Always consume a stored tensor via
    `get_voigt`/`get_tensor`, never by reading `.dev`/`.p` directly.**
  - **Return-mapping kernels stay in Voigt space internally** (`DP.jl`/`J2.jl` read
    `get_voigt(...)` once, run the closed-form/CPA algebra on `SVector`s, wrap once at
    the store) — deliberate, keeps numerics bit-identical to the pre-port code and
    preserves `elast_fast`'s hand-inlined performance.
  - Numerically verified bit-identical predictor math; the storage round-trip itself
    costs ≤1 ulp of the largest component. `test_performance.jl` measured -12.1%
    memory/-5.4% allocs/0.0% time vs. the pre-port baseline (this repo's typed-tensor
    objects allocate less than the bare `SMatrix`/`SVector` fields they replaced). JLD2
    needs no special handling — plain immutables over `T`/`SMatrix`.
  - `get_J2` (the tensor-invariant helper) was deliberately **not** unified with
    `retmap/J2.jl`'s same-named helper (which computes a different thing, `(‖ξ‖,n̂)`) —
    that one was renamed `_yield_normal` instead, so the collision is gone.
- ~~`PointSolidPhase.elast::E`/`AbstractElasticity`/`FiniteElasticity`/
  `LinearElasticity`~~ — **removed**, along with the `E` type parameter. `E` had
  decayed into a pure dispatch tag duplicating what `ST` already encodes once the
  tensor port landed (confirmed zero reads of `mpts.s.elast` before removal). Every
  kernel that dispatched on `E<:FiniteElasticity`/`E<:LinearElasticity` (`elast.jl`,
  `dynamic_relaxation/{fint,update}.jl`) now dispatches on `ST` instead.
- ~~`PointSolidPhase.P::Vector{T2}`~~ — **removed**, zero reads; pressure lives on the
  stress objects themselves (`mpts.s.σᵢⱼ[p].p`, or `get_voigt`'d).
- **`PointSolidPhase.ϵpII` is `Vector{SVector{2,T2}}`**, not `Matrix{T2}` — matches the
  rest of the codebase's per-particle `Vector{SVector}` convention. Since `SVector` is
  immutable, every in-place mutation site (`DP.jl`/`J2.jl`'s CPA updates, `nonlocal.jl`,
  `update.jl`'s bulk resets) reads the whole `SVector`, rebuilds it, writes it back —
  `drucker_prager`'s internal scratch stays a mutable `MVector{2,T}` for the actual
  return-mapping algebra; only the *stored* `mpts.s.ϵpII[p]` is immutable. Measured
  perf-neutral-to-better vs. the old `Matrix` layout (the `+=`-mutation pattern is
  ~25% *faster* as `Vector{SVector}`).
- **`get_voigt` (renamed from `get_vector`) is the single entry point for
  tensor↔Voigt-vector conversion, using the *engineering*-shear convention
  (`γxy=2εxy`)** — required wherever strain interacts with the elastic stiffness
  matrix `Del` (`σ=Del·ε_voigt` only works with the ordinary isotropic stiffness
  matrix under that convention). The old free functions `mutate`/`_mutate` existed
  purely to patch a factor-of-2 mismatch between `get_vector`'s old tensor-shear
  convention and `Del`'s engineering-shear expectation at each call site; both are
  gone now that `get_voigt` uses the correct convention directly.
  `LogarithmicStrain(ϵ::SVector)`/`InfinitesimalStrain(ϵ::SVector)` full-Voigt-vector
  constructors exist too, mirroring `CauchyStress(σ::SVector)`/`KirchhoffStress`;
  `get_voigt(LogarithmicStrain(ϵ)) == ϵ` to machine epsilon. `voigt_of(M::SMatrix)`
  (`elast.jl`) handles the one remaining case with no typed-tensor home (the Jaumann-
  rate correction term `σJ*ω'+σJ'*ω`, not itself "the" particle's stress).
  `retmap/MCRetMap.jl`'s own `mutate(...)` call is dead/unwired code, left calling an
  undefined function if ever exercised — not in scope to fix.
- **Unicode subscripts encode tensor rank — see `conventions.md`.**
  `PointSolidPhase.σᵢ`/`τᵢ` were renamed `σᵢⱼ`/`τᵢⱼ` for exactly this reason (a
  single-index subscript misleadingly read like a vector component on a full 2nd-order
  tensor field). **Gotcha discovered doing this rename**: a blind substring rename also
  caught genuinely-vector *local* variables (`get_voigt(...)` results inside
  `DP.jl`/`J2.jl`/`dynamic_relaxation/fint.jl`), which should have stayed single-index
  — fixed by renaming those specific locals back. Untangling this also surfaced a real,
  separately-introduced bug: `drucker_prager`'s return-mapped stress was briefly being
  written to a variable the function no longer returned, so it silently returned the
  pre-yield stress on every yielding step — verified fixed via a direct unit test (a
  hand-built past-yield stress state now returns an actually-changed result). Lesson:
  a mechanical rename like this needs a semantic check (tensor vs. vector), not just a
  substring match.
- **3D conformity check — done, via `test_workflow.jl` itself.** See
  `.claude/bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md` for the full
  results (149/192 passed overall, 3D not meaningfully less stable than 2D) and the
  related `plast.status` finding.
- **`basis.how`/`basis.ghost` are GIMP-specific concepts living in the generic `basis`
  config section.** `how` is read unconditionally in `update.jl`'s dispatch regardless
  of basis kind; worth moving both onto `GimpBasis` itself (construction-time fields
  co-located with the kind that actually uses them) rather than basis-kind-agnostic
  top-level knobs `bsmpm`/`smpm`/`mlsmpm` silently ignore.
- **`Basis.N`/`∂N` storage rework**: currently plain `Matrix{T2}`/`Array{T2,3}`
  specifically to dodge the StaticArrays allocation-elision limits (see `gotchas.md`)
  — correct and fast, but a workaround rather than a considered data-layout design.
- **`smpm`'s isolated gradient-consistency glitch**: `test_basis.jl` found
  `max|Σᵢ∂Nᵢ|` off by exactly `1.0` at one sweep point for `LinearBasis`, unlike
  `gimpm`'s broadly-degraded boundary PoU. Never root-caused — an error of exactly
  `1.0` smells like a real (if narrow) bug rather than float noise.
- **`DynamicRelaxationSolver`'s fate**: now that `solution` picks between
  `ExplicitSolver`/`ImplicitSolver` only, this third struct is orphaned scaffolding.
  Decide: wire it in as a third `solution` value, or state plainly that it's
  intentionally dormant.
- **Staged path to a poro-hydro-mechanical solution: mechanical solver (exists) →
  standalone fluid-only (hydro) solver → fuse the two.** Validate the hydro solver
  independently before coupling. `PointFluidPhase` (`lagrangian.jl`) needs real fields
  first — currently a literal empty placeholder, `mpts.f` always `nothing`.
- **Thermal solution — fixed and viable.** `thermal_problem` (renamed from
  `ic_thermal`) → `elastoplasm(jld2; workflows=[thermodynamic!])` runs end-to-end
  (verified: temperature stays bounded, cools monotonically from the boundary inward,
  both with `musl=true`/`false`). What was broken: stale `mpts.ϕ∂ϕ` reads across
  several kernels (replaced with the same `basis.N`/`∂Nrow(basis.∂N,...)` caching
  pattern the solid phase uses); a mistyped thermal `p2n!` overload dispatching on
  `mesh::Mesh` instead of `MeshThermalPhase`; no thermal-phase `mapsto` overload
  existed at all (added, mirroring the solid-phase one minus gravity/APIC/finite-
  strain concerns); `mesh.t`/`mpts.t` were unconditionally `nothing`
  (`setup_mesh`/`setup_mpts`/`setup_problem` gained a `thermal::Bool=false` kwarg);
  initial temperature was always zero regardless of config (off-by-one in
  `get_thermal.jl`). Not fixed, left as a known gap: `thermal_problem` defaults to
  `plot.status=false` — the plotting path has no `"T"` field entry and
  `display.jl` carries two incompatible `what_plot_field` dispatch styles, a separate
  piece of follow-up work. Separately found+fixed in the same pass: 3D APIC transfer
  for the *solid* phase had the identical dead-`ϕ∂ϕ` bug plus a mismatched mesh-phase
  type — fixed, but 3D+APIC then hits the pre-existing 3D instability noted in
  `.claude/bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md`; 1D APIC has a
  similar unfixed inconsistency (1D is essentially never exercised in this codebase).
- **Consider a persistent Julia session for iteration (Kaimon/Revise pattern, seen in
  sibling repo FEMTools.jl).** A `slump_problem`→`elastoplasm!` smoke test costs full
  Julia startup+precompile every invocation today (~30-40s), and
  `elastoquasistatic!` is worse (~1 minute) — exactly the kind of loop a persistent
  `Revise`-backed session amortizes. FEMTools.jl's `.agents/`-documented setup
  (`kaimon -r`, exposing `julia_eval`/`julia_list_sessions`/`julia_restart` over MCP,
  restarting only on `Project.toml`/`struct`/include-order changes Revise can't track)
  is a candidate pattern to port if the user sets up the equivalent tooling here — not
  something to install unprompted, since it's an external dependency/process the user
  needs running before a session starts.
- **Add a JLD2 time-series export option, as an alternative to plot-only checkpoints.**
  Every workflow's time loop calls `bake(mpts,mesh,t,solver)` at each checkpoint, but
  `bake` only renders a PNG when `solver.plot.status` is true — there's no way to
  record a field's evolution as raw data today. A natural fix reuses the same
  `checks` cadence and `solver.plot.what`-style field selection to append values into
  a JLD2 dataset at each checkpoint (`export.status`/`export.what`, sitting next to
  `bake()`). Writing the export is only half the feature — a companion loader/plotter
  for the saved time series would need to be built alongside it, not as a separate
  follow-up.
