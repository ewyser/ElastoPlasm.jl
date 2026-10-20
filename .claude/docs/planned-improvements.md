## Planned improvements (not bugs, just known follow-up work)

- **`Problem` — done: one object holding a run.** `Problem` (`types/problem/problem.jl`) bundles
  `mpts`, `mesh`, `basis`, `time`, `solver`, built by `Problem(case)` from a `Case`; it replaced
  `MechanicalProblem` (`mesh`/`mpts`/`time` only, whose docstring kept `Basis`/`Solver` out because
  `basis.which`/`basis.trsfr` vary independently, which does not argue against bundling them for a
  run). The `Case` is deliberately not a field: the built objects already hold all it describes, it
  would drift once `mpts` evolves, and a `Region` closure does not serialise; its TOML beside the
  `.jld2` is the provenance. Kernels keep their `(mpts, mesh, basis, …)` arguments; workflows have
  `Problem` methods and `elastoplasm` falls back to the five-argument form. Possible follow-up: the
  solver mixes configuration (part of the problem) with the kernel table `cairn` (machinery,
  rebuilt from the configuration); rebuilding `cairn` on load would make the stored `Problem` pure
  data. `types/problem/` files load in `boot.jl`'s explicit order, `problem.jl` last (it needs
  `Point`, `Mesh`, `Basis`, `Time`, `Case`).
- **Kernel granularity**: `elast!` is the model to follow — decomposed into small
  functions each doing one specific task, rather than one monolithic kernel body.
  `update.jl`'s `elastoplast`/`elasto` dispatch functions (branch on
  `basis.which`/`stab.locking`/`nonloc.status` inline) would benefit from the same split.
- **Typed constitutive-model abstraction — done.** `mpts.s.cmp::Vector{CM} where
  CM<:AbstractConstitutiveModel` (`constitutive.jl`) bundles the *static* elastic+
  plastic material constants (`Gc`, `Kc`, `Hp`, `c₀`, `cᵣ`, `ϕ₀`) into one typed
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
- **Typed strain/stress tensor storage — done, as a scalar plus a 3D Voigt deviator.** `mpts.s.ϵᵢ`/
  `.ϵn` (`LogarithmicStrain{T}`/`InfinitesimalStrain{T}`, picked by `material.elastic`) hold `vol::T`
  (the trace) and `dev::SVector{6,T}`; `.σᵢ`/`.σn`/`.τᵢ` (`CauchyStress{T}`/`KirchhoffStress{T}`)
  hold `p::T` (`-tr/3`, positive in compression) and `dev::SVector{6,T}` (`tensor.jl`, loaded
  before `lagrangian.jl` in `boot.jl`'s explicit list; the elastic step, i.e. the Hencky push-forward,
  the laws and their inverse, sits in `elasticity.jl` right after it). `dev` is trace-free Voigt
  `(xx, yy, zz, yz, xz, xy)` in every dimension: a 2D problem is plane strain (`yz = xz = 0`, `zz`
  carried); 1D is uniaxial strain. Strain shears are **engineering** (`γ = 2ϵ`).
  - **All constitutive code runs on the split.** The Hencky laws: `p = -Kc·vol`,
    `dev_τ = 2Gc·dev_ϵ`, the engineering shears halved to tensor components (`_shear_weights`); their inverse `_elastic_strain`
    rebuilds the strain on yield in both return maps; the hypoelastic increment is the same law
    applied to `ϵ`. Drucker–Prager reads and returns `(p, dev)`. The von Mises CPA moves only
    `dev` (`Del·n̂ = 2G·n̂` for a deviatoric `n̂`, so the pressure is unchanged).
  - **Kinematics meets the split in one place each way:** the 3×3 constructors (`vol = tr`,
    `dev = M − tr/3·I`; the push-forward, the hypoelastic strain increment, the Jaumann term) and
    `get_tensor` (`dev ± …·I`; P2G, DR `fint`). There is no Voigt-vector constructor and no
    projection in the constructors: a `(p, dev)` built directly must already be trace-free, which the
    laws and return maps are by construction.
  - Other readers: `get_tensor(t, Val(D))`, `get_J2`/`get_τII`, `eigen(ϵ)`; pressure and
    volumetric strain are read as `.p`/`.vol`. `get_voigt`, `pressure` and `volumetric` were
    removed once nothing needed a full Voigt vector (1D P2G and `test_column` read
    `get_tensor(σ)[i,i]`, bit-identical).
  - Split storage replaced single `v::SVector{6}` storage (the user preferred the separate
    pressure/deviator concerns). Reproduced it to round-off: DR blocks ≤3e-13 m with identical
    iteration counts, explicit slumps ≤5e-9 m, hypoelastic slumps ≤7e-13 m, column errors to 1e-12,
    `test_column` 4/4. Cost, best of 5 alternating runs: `elast!` 59.6 → 60.3 µs, VM `retmap!`
    490 → 480 µs, DP `retmap!` ~143 → 145 µs. Moving the return maps and the conversions onto the
    split (no `Del` solve, no Voigt round trip) again reproduced results to round-off, and cut the
    cost (median of 5 alternating runs): `elast!` 59.7 → 55.9 µs, VM `retmap!` 480 → 161 µs, DP
    `retmap!` 146 → 40 µs.
  - **3D checked against `dev`** (slump `[64.16, 16.04, 16.04]` m, `[40,10,10]` elements, 11,820
    particles): Drucker–Prager matches `dev` to round-off (positions ≤3e-13 m, stress ≤2.4e-13
    relative). Von Mises differs by up to 4.5 m, entirely from the radial-return fix `fb794f1`: at
    that commit 3D von Mises already matches this branch to 5e-13 m, so no later storage change
    alters 3D results.
  - Constitutive code runs once, in 3D: `elast!` keeps the kinematics D×D (the hypoelastic strain
    increment carries zeros, the Jaumann-rotated stress carries `σzz`, both via the carry
    constructor `S(M, t)`);
    both return maps work on the 3D split; Drucker–Prager keys its
    2D fit on the model's `D`. Hencky's push-forward (`_trial_elastic_strain`) works on the
    in-plane D×D block only and carries the out-of-plane components through unchanged: with
    `ΔF_zz = 1` and no out-of-plane shear, `b_zz = exp(2ϵzz)`, so `ϵzz` is exactly unchanged.
    This replaced a 3×3 push-forward (block-split `eigen`, 3×3 matrix products): `elast!` 67 →
    60 µs (pre-refactor 56 µs; best of 5 alternating runs), results equal to round-off.
  - History: the 3D storage replaced a `(p|vol, dev::SMatrix{D,D})` split, which in 2D could not hold `zz`
    (so plane strain needed a separate `mpts.s.ϵzz` field and 2D-only return methods) and had
    constructors disagreeing on `tr/2` vs `tr/3` (so its `dev` was not trace-free). The switch (branch
    `feat-refactor-tensor-storage`) reproduced the previous results to round-off — DR to 1e-13
    with identical iteration counts, explicit slumps to 7e-9 m, `test_column` errors to 6e-13,
    `test_workflow` matrix 145/192 before and after with the same failing set — at net −336
    lines. Intended change: hypoelastic plasticity now uses full plane-strain invariants
    (plastic slumps move up to 0.77 m). Cost: `elast!` 56 → 74 µs on a 591-particle slump
    (fresh state), since reduced to 60 µs by the in-plane push-forward; `elast!` is a few % of a time step. `get_P` became the true mean pressure
    (`σzz` included) instead of the old "mean of stored Voigt components". `.jld2` files written
    before the switch don't load.
  - `get_J2` (the invariant) and `retmap/J2.jl`'s `_yield_normal` (`(‖ξ‖, n̂)`) are different
    helpers on purpose.
  - **`Del` is not stored; it is built on the fly.** `get_Del(cmp) = Kc·vol + 2Gc·dev`
    (`elasticity.jl`, `StressStrainStiffness{T}`: the 6×6 volumetric/deviatoric Voigt projections,
    `½` on the shear diagonal for engineering strain). `DruckerPrager`/`VonMises` are now
    `{T2,D}` with no stored `Del`, and `get_elastic_stiffness` became `get_elastic_moduli(E, ν) ->
    (Kc, Gc)`. Matches the old stored `Del` to ≤2.5e-16 relative; results: Drucker–Prager
    bit-identical, von Mises to round-off (≤1e-13 m DR with identical iterations, ≤8e-12 m explicit);
    per-particle `cmp` 328 → 40 bytes; von Mises `retmap!` time unchanged (~500 µs, within noise).
    A future non-constant stiffness (e.g. improved Hencky's tangent bulk modulus) only needs a
    different `Kc` passed to the same projections. Since the constitutive code moved onto the split, nothing calls
    `get_Del`; it is kept on purpose for future use (e.g. a consistent tangent).
- ~~`PointSolidPhase.elast::E`/`AbstractElasticity`/`FiniteElasticity`/
  `LinearElasticity`~~ — **removed**, along with the `E` type parameter. `E` had
  decayed into a pure dispatch tag duplicating what `ST` already encodes once the
  tensor port landed (confirmed zero reads of `mpts.s.elast` before removal). Every
  kernel that dispatched on `E<:FiniteElasticity`/`E<:LinearElasticity` (`elast.jl`,
  `dynamic_relaxation/{fint,update}.jl`) now dispatches on `ST` instead.
- ~~`PointSolidPhase.P::Vector{T2}`~~ — **removed**, zero reads; pressure comes from the
  stress objects themselves (`mpts.s.σᵢ[p].p`).
- **`PointSolidPhase.ϵpII` is `Vector{SVector{2,T2}}`**, not `Matrix{T2}` — matches the
  rest of the codebase's per-particle `Vector{SVector}` convention. Since `SVector` is
  immutable, every in-place mutation site (`DP.jl`/`J2.jl`'s CPA updates, `nonlocal.jl`,
  `update.jl`'s bulk resets) reads the whole `SVector`, rebuilds it, writes it back —
  `drucker_prager`'s internal scratch stays a mutable `MVector{2,T}` for the actual
  return-mapping algebra; only the *stored* `mpts.s.ϵpII[p]` is immutable. Measured
  perf-neutral-to-better vs. the old `Matrix` layout (the `+=`-mutation pattern is
  ~25% *faster* as `Vector{SVector}`).
- **Strain deviators use the engineering-shear convention (`γxy=2εxy`)**, what `Del` acts on. The
  old `mutate`/`_mutate` patches for a factor-of-2 mismatch, `voigt_of` and `get_voigt` (formerly
  `get_vector`) are gone.
- **Unicode subscripts reflect storage — see `conventions.md`.** Stress/strain fields are
  `σᵢ`/`τᵢ`/`ϵᵢ` since the switch to Voigt storage. (History: they were once renamed `σᵢⱼ`/`τᵢⱼ`
  under an earlier rule tying subscripts to physical tensor rank.) **Gotcha discovered during
  that earlier rename**: a blind substring rename also
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
- **`basis.how` and the alternative GIMP domain updates — removed.** GIMP's particle domain
  is always updated by `Uᵢᵢ` (stretch `U` from `F`, `domain.jl`). The six alternatives
  (`undeformed`, `detFᵢᵢ`, `Fᵢᵢ`, `detΔFᵢᵢ`, `ΔFᵢᵢ`, `ΔUᵢᵢ`) and the `basis.how` knob meant to
  pick between them were only ever referenced from a commented-out block, so `how` was read
  nowhere. `test_workflow.jl`'s GIMP case passed `how="undeformed"` while actually running
  `Uᵢᵢ`. (`basis.ghost` no longer exists either.)
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
  a finely meshed `elastoquasistatic!` run adds more — exactly the kind of loop a persistent
  `Revise`-backed session amortizes. FEMTools.jl's `.agents/`-documented setup
  (`kaimon -r`, exposing `julia_eval`/`julia_list_sessions`/`julia_restart` over MCP,
  restarting only on `Project.toml`/`struct`/include-order changes Revise can't track)
  is a candidate pattern to port if the user sets up the equivalent tooling here — not
  something to install unprompted, since it's an external dependency/process the user
  needs running before a session starts.
- **Time-series output — done, as HDF5** (`output` solver section, `src/home/api/io/output.jl`,
  `src/home/plot/output_plots.jl`; see "Recording a run" in `operations.md`). Follow-ups:
  whole-model series (kinetic energy, run-out, max displacement; one value per frame, `/series/…`);
  an XDMF sidecar so ParaView opens the `.h5` as an animated point cloud; fields of the fluid phase
  once it has some; `plot_output` in 3D shows the x–z projection of all points (a slab or slice option
  would help); the figures share the live plots' cramped sizing.
- **Decide when `dump/` is emptied.** It happens only when the package precompiles (see
  `conventions.md`), not at each `using` as the docs used to say: either move `rootflush` to
  `__init__` (every load, the documented behaviour) or drop it (outputs then accumulate until removed
  by hand). The user's call.
- **Improved Hencky's porosity is its own elastic-strain measure — decided, keep.**
  `_trial_elastic_stress_improved` (`elasticity.jl`) and `_Ktan` (`explicit/get.jl`) compute
  `n = 1-(1-n₀)/exp(ϵᵥᵉ)`, i.e. Eqs. 22-23 of Pretti et al. (2024) written in the stored
  elastic log strain, rather than reading `mpts.n`. This keeps `n` and `ϵᵥ` on the same
  measure, which is what makes Eq. 34 the derivative of Eq. 33. `mpts.n`/`ρ`/`Ω` come from the
  raw `J = det(F)` in `deform.jl`, so the two porosities differ whenever F-bar is on (the
  default — `ΔJp!` rescales `ΔFᵢⱼ`, never `Fᵢⱼ`/`J`) or plastic flow is dilatant. The paper
  has no F-bar, so the combination is off-paper. Follow-up: the `elasticity.jl` docstring still
  says the two agree "while plastic flow is isochoric", which leaves out F-bar. Related, open:
  `deform.jl` computes `ρ = (1-n)ρ₀` from the *unclamped* `n` but stores `mpts.n` clamped to
  `[0,1]`, so once compaction crosses Pretti's Eq. 30a (`J < 1-n₀`), `mpts.n` sits at 0 while
  `ρ` keeps growing past `ρ₀`. `ρ·Ω` (the solid mass) stays exact either way.
- **DP apex `ϵpII` increment misses the shear-flow part** — deferred follow-up of the fixed
  apex return (`.claude/bug/fixed/dp-apex-return-wrong-pressure.md`). Needs a corner-return
  source (de Souza Neto, Perić & Owen 2008 §8.3) added to `refs/` first.
- **Dead or drifting code worth removing, one small verified step at a time:** no `#= … =#` blocks
  are left under `src/`, and `MeshSolidPhase.oobp`/`.p` (the removed u–P path) are removed. Still
  open: `PerfectlyElastic` and `DynamicRelaxationSolver` above.
- **`_fast` kernels and the `perf` config section — removed.** `deform_fast`/`elast_fast`
  (hand-unrolled 2D/3D scalar versions, a legacy port of the original C code) were only
  selected by `perf.status`, which also forced `material.elastic="hypoelastic"` and
  `nonloc.status=false`. Measured before removal: 1.8–2.7× faster per kernel, but
  `deform!`+`elast!` are only ~3–5% of a timestep, so the whole-step cost of removal is small;
  `deform_fast` had also drifted into a wrong porosity update (`n = 1-(1-n_prev)/J` with the
  *total* `J`, no clamp). Default-path results were bit-identical before/after removal.
  Possible follow-up: speed up the generic `deform` (its `MMatrix` accumulation is the likely
  gap), measured, rather than reintroducing a second implementation.
  The removal also silently dropped the hypoelastic Jaumann rotation, which only the `_fast`
  kernels applied (the generic term was identically zero); see
  `.claude/bug/fixed/hypoelastic-jaumann-rotation-is-zero.md`.
- **GPU execution is unwired — decide whether it's a goal.** `get_solver` resolves
  `backend.exec` from `backend.select`, but nothing reads it: every kernel factory in
  `init_ignite`/`init_mapsto`/`init_update` and every `sync(...)` call is a literal `CPU()`
  (61 occurrences under `explicit/` and `common/`). `get_dt` also reduces over particles in a
  host loop with scalar indexing, which device arrays don't allow. Wiring it is broad,
  mechanical work (thread one backend object through the factories and syncs, verified
  bit-identical on `CPU()` first), only worth doing if GPU runs are actually wanted.
- **`@adapt_struct PointSolidPhase` likely can't rebuild the struct on GPU, unverified.**
  `Adapt` reconstructs a struct from its adapted field values, which needs every type parameter
  to be inferable from the fields. `PointSolidPhase`'s `T1` and `EL` appear in no field.
  (`Point` had the same problem through the phantom `TS` parameter until the type-parameter
  cleanup removed it.) Not testable on the CPU-only setup used so far.
- **`ExplicitSolver`/`ImplicitSolver` are two identical 15-parameter structs kept in lockstep
  by hand**, one parameter per config section. A single solver struct holding the config
  `NamedTuple` as one field plus an explicit/implicit tag would remove most of that, at the
  cost of touching every `solver.xxx` access.
- **`test_workflow.jl` doesn't cover `material.elastic="improved_hencky"`** — add it once
  `.claude/bug/known/improved-hencky-explicit-admissibility-crossing.md` is fixed.
- **Lighter `src/home/init/` — done** (branches `refactor-lighter-init` → `fix-grf-point-evaluation`
  → `refactor-init-problem-spec`). Every problem is a `Case` built by `setup_simulation`; the
  six hand-written pipelines and the `get_*` generators are gone (see "Describing a case" in
  `operations.md`). Material constants live in `MATERIAL` (`setup_material_constants` removed, `n₀`
  a material constant), `PointData` replaced the implicit generator NamedTuple and
  `build_solid_phase`'s 11 positional arguments, `get_solver` merges per field, a 3D column works.
  Verified bit-identical setup by setup with `test_setup.jl` (a stored digest reference, which
  replaced the scratch harness used for the earlier steps). The `Geometry` type is gone too:
  `setup_mesh(L, nel, solver)` builds the node grid once, for any dimension (`mesh_coords`, which
  errors if a floating-point range misses a node), instead of a `Geometry` whose counts `get_coords`
  recomputed and `MeshProperties` stored again. Earlier steps: one `mpts_populate` grid
  construction, the slump slope as one half-plane test.
  `init/` docstrings describe behaviour (Step 8). Remaining:
  - **Step 5b, per-property random fields (deferred until a second random property is needed):**
    `case_points` would take a per-property `grf.fields` config and call `property_field` for each
    listed property. Open questions: truncation at a floor vs lognormal mapping for positive-only
    properties (needs a reference in `refs/`, e.g. Fenton & Griffiths 2008); whether zone values
    override or centre the random value (today a zone setting `c0` overrides it); independent vs
    cross-correlated fields (shared harmonics).

- **Setup cost — `e2e` fixed; IC plot and export are what's left.** `e2e` (`init/basis/e2e.jl`)
  filled a `SparseMatrixCSC` entry by entry; in 1D/3D each row insert shifted all later storage, so
  the build was quadratic (3D slump `[40,10,10]`: 2.4 s; `[60,15,15]`: 28 s; 32,000 elements:
  174 s). It is now one `sparse(I, J, V)` call over `CartesianIndices`, one code path for every
  dimension: 4 ms and 17 ms for those meshes, near-linear up to ~1e6 elements. Same sparsity
  pattern (the only thing `nonlocal.jl` reads); 2D and 3D slumps bit-identical. Whole setup
  (compilation excluded): 3D `[40,10,10]` 4.6 → 2.7 s, `[60,15,15]` 36.6 → 8.5 s. What remains is
  the IC plot (`slump_problem` with `plot.status`, 2.1 s / 5.2 s) and the `.jld2` export (0.3 s /
  2.1 s); both linear, not yet looked at.
- **Random field (`GRF.jl`) open points**, none blocking: (a) the per-axis extension of the
  Gaussian covariance is derived in-repo, not stated in Räss et al. (2019); Sabelfeld (1991) is
  a candidate source. (b) Point sampling at particle centres vs local averaging over the
  particle volume (Vanmarcke; Fenton & Griffiths), which matters more for the rougher
  exponential covariance and when particle spacing approaches the correlation length.
- **Mesh variables can't be plotted through `plot.what`.** `get_mesh_variable_config()`'s only entry
  (`"m"`) has `data=get_epII` and no `what` field, while `what_plot_field(mesh, opts)` branches on
  `opts.what`, so a `(; mesh=...)` entry from `get_variable_plot_options()` would throw.

- **Usability findings (second `agent-user` pass), not yet addressed:** a mistyped nested config
  key (e.g. `material=(;elastik=...)`) replaces the whole section and fails in `get_solver` with
  a cryptic "NamedTuple has no field" error; reading typed tensors (`.p`/`.dev`, `get_tensor`) and the
  `dump/` flush on `using` are undocumented in `README.md`/`docs/src`; `worflow.jl` is misspelled
  in both solver folders.

- **Implicit dynamic-relaxation MPM (Sutcliffe et al. 2026) — Algorithm 1 and plasticity
  done**, see `architecture.md`. Remaining, in order:
  1. Small cuts — addressed for DR by fixing low-solid-fraction nodes (see "Nodal solid fraction"
     below). Original notes: the paper's aggregation (Alg. 2) is **not** planned — it only applies to
     linear-element bases, not `bsmpm`. The intended route is to evolve a solid fraction and
     drop empty or near-empty nodes from the residual assembly, after
     `refs/ElastoPlasm/fractional_step/s00366-026-02369-3.pdf`. Symptom it should remove:
     occasional DR iteration spikes when the column's top surface crosses a cell edge, and
     the collapse's 10k–1M-iteration load steps.
  2. Open reading choices, recorded when implemented: the starter `v₋½ = −½M̃⁻¹fᵉˣᵗₙ` is
     taken literally (full load, not the increment); `ε_III` in 2D includes the zero
     out-of-plane strain (the 3D strain's `zz` component); `‖R‖` is over free DOFs only; the
     paper's `ρ_y` is read as a Lode radius (`√(2J₂)`), which is what reproduces Fig 17a.
- **2D plane-strain invariants — done for every constitutive path.** Stress and strain are full
  3D Voigt vectors (see "Typed strain/stress tensor storage"), so von Mises and Drucker–Prager,
  finite and small strain, include `σzz` in J₂ and the mean stress (Huang et al. 2015, Eqs 15-16).
  Verified: von Mises return radial and matching a 3D return, K₀ yield onset between 60 and
  61 kPa (`κ = √2·20 kPa`, `ν = 0.3`), DP on the yield line with pressure unchanged (`ψ = 0`) and
  apex return under tension. `σzz` follows the elastic law of the full 3D strain.
  - **DP's 2D parameter fit is not yet grounded in `refs/`.** 2D keeps the plane-strain match
    to Mohr–Coulomb (`q_ϕ = 3 tan ϕ/√(9+12 tan²ϕ)`, `k_ϕ = 3c/√(9+12 tan²ϕ)`, FLAC theory
    manual), chosen because it is meant for full invariants under plane strain; Huang et al.
    (2015) only give the inner-cone fit used in 3D. Add the FLAC manual (or an equivalent
    source) to `refs/` and check the fit's assumptions (it is derived for associated flow; the
    code runs `ψ = 0`).
- **Nodal solid fraction, fixing small-cut nodes in DR — done, DR only.** `src/home/core/common/sfrac.jl`,
  after Xie, Navas & López-Querol (2026, *Eng. Comput.* 42:146, `refs/ElastoPlasm/fractional_step/
  s00366-026-02369-3.pdf`): nodal volume `V_I = Σ S_I,gp V_gp` over a completely filled grid (Eq 37,
  quadrature points laid out like the material points, the problem's own basis; interior nodes get
  exactly `hᴰ`, `Σ V_I` = domain volume), and relative solid fraction
  `r_I = Σ_p S_Ip V_p (1−n_p)/(1−n₀_p) / V_I` (Eq 36 normalised by the initial solid fraction, so 1
  inside the material; it equals `Σ S_Ip V⁰_p / V_I`). `relax!` fixes (`Δu = 0`) nodes with
  `0 < r_I < stab.sfrac` (default `0.003`, `0` = off) for the load step — the paper's Eq 45
  criterion, whose action there is `p = 0` on the water free surface; for a displacement solve the
  analogue is fixing, which the user chose.
  - **Effect** (§4.3 VM collapse, Lode reading, 2² MPs/cell, 10 load steps): 1,204,516 DR iterations
    (one load step capped at 1e6) → 16,058; the √2-reading VM block 69,074 → 10,199; the DP block
    41,384 → 6,284. `test_column` unchanged (4/4).
  - **Bias: fixing holds the free surface back, in proportion to the threshold.** Run-out vs `ε`:
    0 → 5.74 m, 0.001 → 5.74, 0.003 → 5.69, 0.01 → 5.46, 0.03 → 4.13 (iterations 1.2M, 49k, 16k,
    15k, 11k); `0.003` is the chosen trade-off (<1% run-out bias). Locally the bias can be larger:
    in the √2-reading VM block one particle's stress changed by 42% of the peak stress (positions by
    ≤3.4 cm), most likely a surface particle beside a fixed node.
  - **Not used by the explicit solver**: tried and reverted. Fixing nodes in a flowing explicit run
    made the slump crash after ~3,400 steps (a particle's `b` lost positive-definiteness) —
    hypothesis, not proven: surface particles get a spurious velocity gradient from the fixed nodes
    every step, which accumulates. The explicit runs were stable without it.
  - Follow-ups: an extrapolated (not zero) displacement at fixed nodes would remove the surface bias
    but needs a reference; this stays open alongside the paper's aggregation, which is not planned.
- **`test_performance.jl` rewrites its tracked baseline on any machine not in it.** The baseline
  (`test/dataset/performance_baseline.jld2`) holds one AMD Ryzen entry; elsewhere the test saves the
  current run as that machine's entry (a change to a tracked file) and then compares the run with
  itself (0.0% everywhere). Restore the file after running it on another machine, or keep machine
  baselines out of git.
- **`test_performance.jl` benchmarks `elast!` on NaN state.** The suite runs the mapping and
  deformation kernels repeatedly on an unsettled problem without resetting it, so `elast!` sees
  NaN input, whose 3×3 symmetric `eigen` is ~4000× slower than a valid one. Its `elast` number is
  therefore meaningless since the 3D storage; benchmark each kernel on a valid, reset state (as
  the scratch `elastbench.jl` did: 56 → 74 µs), then regenerate the baseline deliberately.
