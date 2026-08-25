# ElastoPlasm.jl — Notes for Claude

Working notes on this repo's architecture and conventions, kept up to date as the
codebase evolves. Read this first when picking up work here.

## Project shape

MPM/geomechanics solver in Julia. Files are loaded via a custom recursive includer
(`src/boot/include.jl`, `collect_and_include_jls`) that walks `src/` and includes every
`.jl` file — files in a directory load before its subdirectories, both sorted
alphabetically, so load order is deterministic. A failed include only `@warn`s, it
doesn't crash — so a stale/broken file can silently go unused; always check for
`@warn "Failed to include..."` after `using ElastoPlasm` when touching files here.

## Core types

- `Mesh{T1,T2,D}` — Eulerian background grid. No longer carries `NN` or connectivity
  fields (`e2n`/`e2e` moved out, see `Basis` below).
- `Point{T1,T2,D,CM<:AbstractConstitutiveModel,TM,TV,TS,ST<:AbstractStrain,SC<:AbstractStress,SK<:AbstractStress}`
  — material points (Lagrangian). The trailing `ST`/`SC`/`SK` are the typed
  strain/stress storage on `mpts.s` (see "Typed strain/stress tensor storage" under
  Planned improvements); `mpts.s.σᵢ[p]` returns a `CauchyStress`, **not** an
  `SVector` — read it with `get_voigt(...)`. No longer carries `NN`, `B`, or connectivity
  (`p2n`/`p2e`/`e2p`/`p2p` moved out). `mpts.x :: Vector{SVector{D,T2}}` — NOT a
  matrix; index fields with `getindex.(x, i)` or iterate, never `x[i, :]`.
  `mpts.s.cmp::Vector{CM}` is the per-particle constitutive-model bundle (see "Typed
  constitutive-model abstraction" under Planned improvements). The `E<:AbstractElasticity`
  type param (and `AbstractElasticity`/`FiniteElasticity`/`LinearElasticity` themselves)
  were removed entirely — `E` had decayed into a pure dispatch tag duplicating what `ST`
  (`LogarithmicStrain` vs `InfinitesimalStrain`) already encodes once the tensor port
  landed, so every kernel that used to dispatch on `E<:FiniteElasticity`/
  `E<:LinearElasticity` (`elast.jl`'s `elast`/`elast_fast`, `dynamic_relaxation/
  {fint,update}.jl`) now dispatches on `ST` instead. The parallel `R<:AbstractRheology`
  type param/`rheo` field (a deprecated earlier attempt at the same per-particle
  constitutive-data problem `cmp` now solves) was retired entirely — see "Typed
  constitutive-model abstraction" below.
- `MechanicalProblem{T1,T2,D,CM,TM,TV,TS,ST,SC,SK} <: AbstractProblem{T1,T2,D,SP}`
  (`src/boot/needs/types/concrete/problem/problem.jl`) — bundles `mesh::Mesh`+
  `mpts::Point`+`time::Time` as the IC-defining part of a simulation, built via
  `setup_problem` (see "`Problem` type decoupling..." under Planned improvements for
  the full design and what's still a thin-wrapper simplification). Deliberately
  excludes `Basis`/`Solver` — those stay independent types, not part of `Problem`.
- `Basis{T1,T2,D,NN,K<:AbstractBasis}` (`src/boot/needs/types/concrete/basis/basis.jl`) —
  independent struct owning ALL mesh/point connectivity: `e2n`, `e2e`, `p2n`, `p2e`,
  `e2p`, `p2p`, plus `kind::K` (dispatches shape-function evaluation) and `NN` (nodes
  per element, only known once you fix mesh + basis kind, hence its own type param
  rather than living on `Mesh`/`Point`). Also owns `type` (per-axis node boundary-layer
  classification, `BSplineBasis`-specific) and the shape-function cache `N`/`∂N` (see
  below). `T2` was added late — every field is otherwise `T1`-typed, so adding it was a
  mechanical sweep of every `Basis{T1,...}` call site (~25 files); it exists so `N`/`∂N`
  can live on `Basis` at all.
- Construction order is always **mesh → material points → basis** (basis depends on
  both). Function argument order is always **`(mpts, mesh, basis, ...)`** — this
  ordering is deliberate: it signals that the mpts↔mesh relationship is governed by
  basis.
- Shape-function evaluation is `eval_basis(mpts, mesh, basis, ip, nn)` (renamed from the
  old bare `basis(...)` to avoid colliding with the `Basis` type/struct field name).
  Concrete kinds live in `basis/{bsmpm,gimpm,smpm,mlsmpm}.jl`, each defining its own
  `stencils::NTuple{D,UnitRange{T}}` field directly on the kind struct (no separate
  `stencil_range`/`Neighbs` indirection).
- **Shape values/gradients are cached, not recomputed per call site.** `Basis.N` /
  `Basis.∂N` (plain `Matrix{T2}` / `Array{T2,3}`, sized `(NN,nmp)` / `(NN,D,nmp)` — see
  below for why not `Vector{SVector}`/`Vector{SMatrix}`) are populated once per particle
  per timestep by the `shpfun!` ignite kernel (`src/home/core/common/shpfun.jl`), which
  loops `nn ∈ 1:NN` calling the kind-specific `eval_basis` once each. Every P2G/G2P/
  update kernel then just reads `basis.N[nn,p]` / `∂Nrow(basis.∂N,nn,p,Val(D))` (the
  latter a small helper in `basis.jl` — extracts a gradient row via a `Val`-unrolled
  `ntuple`, since `M[nn,:]` on an `SMatrix` allocates for a runtime `nn`) instead of
  calling `eval_basis` inline. `dynamic_relaxation`'s ~14 call sites are the one
  exception — that path has no ignite phase, so they still call `eval_basis` directly
  per `(ip,nn)`.
- `MLSBasis` (`basis/mlsmpm.jl`) is a fourth kind implementing structured-grid Moving
  Least Squares shape functions (Cao et al. 2025, *Comput. Mech.* 75:655–678, §2.3.1 —
  only the structured-grid case; the paper's unstructured-tessellation machinery does
  not apply here). Unlike the other three kinds it needs no boundary node-type
  correction (`basis.type`) — a per-particle moment-matrix projection restores
  partition-of-unity/consistency everywhere on its own — but it does need a
  once-per-particle `(D+1)×(D+1)` matrix build+inversion, which is why the `shpfun!`
  caching above exists: without it, MLS's `eval_basis` would rebuild that matrix from
  scratch at every one of the ~14 call sites instead of once. `basis.which` and
  `transfer.trsfr` (P2G/G2P transfer scheme: std/tpic/apic) are fully orthogonal — any
  basis kind pairs with any transfer scheme, including `mlsmpm`+`std` (MLS is usually
  *presented* alongside APIC in the literature, but nothing here couples them).

## Explicit vs dynamic_relaxation solvers

- `src/home/core/common/` — kernels/wrappers shared by **both** solver paths:
  `tplgy.jl` (`p2e2n`, topology), `shpfun.jl` (`shpfun!`/`Dij_nd`), `ignite.jl`
  (`init_ignite`/`ignite()`). These used to live under `solver/explicit/ignite/`, which
  was misleading — both `elastodynamic!`/`elastoplastic!` *and* `elastoquasistatic!`
  call the same `ignite(mpts,mesh,basis,solver::AbstractSolver)` wrapper, dispatching
  on the *abstract* solver type since `dynamic_relaxation` doesn't require any
  particular concrete `solution` value (unlike `elastodynamic!`/`elastoplastic!`,
  which are typed to `ExplicitSolver` specifically — see `solution` in "Controlling
  solver behaviour" below).
- `src/home/core/solver/explicit/` — the primary, actively-used path
  (`elastodynamic!`, `elastoplastic!`). Kernel dispatch table lives on
  `solver.cairn` (a `NamedTuple`, dot-accessed: `solver.cairn.ignite.tplgy!`,
  `solver.cairn.mapsto.map.p2n!`, `solver.cairn.update.deform!`, ...), built by
  `init_ignite`/`init_mapsto`/`init_update`/`init_implicit` in `get_solver.jl`. Every
  kernel/function in this path takes `basis::Basis` as its 3rd positional arg after
  `(mpts, mesh, ...)`.
- `src/home/core/solver/dynamic_relaxation/` — quasi-static path (`elastoquasistatic!`,
  plus the not-fully-wired u-P variant `elastoquasistaticuP!`/`elastouP!`). Converted to
  the same `Basis`-threading convention. Runs end-to-end through the same
  `slump_problem`/`elastoplasm(jld2; workflows=[...])` pipeline as the explicit solver — but
  is much slower to smoke-test (a full `elastoquasistatic!` run took ~5–8 minutes in
  practice, vs. seconds for `elastodynamic!`/`elastoplastic!`), so budget for that when
  verifying a change touches this path.

## Known bugs

Bugs discovered but not yet fixed, kept here so they don't get rediscovered from scratch.
**When picking one up: create a new branch off the current branch, named so the fix is
identifiable (e.g. `fix-volumetric-locking-zero-mass-guard`), rather than fixing in place
on an unrelated branch.**

- ~~`basis.which` ∈ {`"gimpm"`, `"mlsmpm"`} combined with `stab.locking=true` crashes
  with `InexactError: trunc(Int64, NaN)`~~ — **fixed** (commit `455fc35`, "fix: Add
  `if...end` check for zero nodal mass"): `volumetric.jl`'s `ΔJp` now has the same
  `iszero(mesh.s.m[no])` guard the velocity kernels already had. Kept here only so the
  fix isn't rediscovered as a new bug.
- ~~`mpts.s.Δλ` is never written under `plast.constitutive="J2"`, silently disabling
  non-local regularization~~ — **fixed** (branch `feat-typed-strain-stress-tensors`).
  `J2.jl`'s `finite_J2`/`infinitesimal_J2` computed `Δλ` as a purely loop-local
  quantity feeding `Δσ`/`γ0`, and never assigned `mpts.s.Δλ[p]` at all. Since
  `setup_mpts` zero-initializes it, it stayed `0.0` forever — and `nonlocal.jl` gates
  *every* one of its three branches on `mpts.s.Δλ[p] != T2(0.0)`, so non-local
  regularization never activated under J2 no matter what `solver.nonloc.status` said.
  Fix: accumulate `Δλ` across the CPA iterations, store the sum after the loop, and
  zero it in the elastic (`f ≤ 0`) branch. Measured on a 591-particle `slump_problem`
  J2 + `nonloc.status=true` run, before → after: `Δλ != 0` on **0 → 424** particles
  (finite) and **0 → 470** (infinitesimal); `ϵpII[2,:] != 0` (the non-local field)
  **0 → 460** and **0 → 474**; `max ϵpII[2]` `0.0 → 3.43e-2` / `6.72e-2`. DP runs were
  never affected. Kept here so the fix isn't rediscovered as a new bug.
- ~~`finite_DP` never resets `mpts.s.Δλ[p]` on an elastic step~~ — **fixed** (same
  branch). `infinitesimal_DP` zeroes `Δλ` unconditionally before the yield check and
  then conditionally overwrites it; `finite_DP` only ever wrote it *inside* its
  `if Δλ > 0` guard, so a particle that yielded once stayed permanently flagged
  "yielding" for `nonlocal.jl`'s gate. (The zeroing had clearly been intended: it was
  present in a large commented-out dead block in that kernel, which has now been
  deleted along with the fix.) Same one-line `mpts.s.Δλ[p] = T2(0.0)` as its sibling.
- ~~`_kirchoff_stress`'s 3D Voigt shear ordering disagreed with every consumer~~ —
  **fixed** as a side effect of the tensor port (same branch). The old
  `_kirchoff_stress(::SMatrix{3,3},...)` returned `(τxx,τyy,τzz,τxy,τxz,τyz)`, i.e.
  slots 4/5/6 = xy/xz/yz — but `std_p2n`/`tpic_p2n`/`apic_p2n`'s 3D kernels, and
  `mutate(...,:voigt)`, all use 4/5/6 = yz/xz/xy. So 3D finite-strain runs assembled
  internal forces from transposed shear components. `get_voigt` uses the consumer
  convention, so this is now consistent; note it means **3D finite-strain results change**
  (2D is unaffected — it has only one shear slot). 3D was already flagged as untested
  territory (see "3D conformity check"), and this was not verified against a reference
  solution, only made self-consistent.
- **`dynamic_relaxation`'s u-P path throws if exercised.** Probable cause:
  `implicit.jl`'s `pt_solve_uP!` calls `instr.cairn.implicit.fint_p2n!(...)`, a key
  `init_implicit` never registers. Out of scope until someone designs what `fint_p2n!`
  should do; `elastoquasistatic!` itself doesn't reach this code path, so it isn't
  blocked by this.
- **`dynamic_relaxation`'s plain (non-u-P) path only actually supports
  `strain.deform="finite"`, despite a comment in `implicit.jl` claiming the opposite.**
  `pt_solve!`/`relax` always call `oobf_assembly!(mpts,mesh,basis,g,dt,Kc,Gc)` — 4
  trailing args, the shape `fint.jl`'s `ST<:LogarithmicStrain` (finite-strain) method
  expects. Under `strain.deform="infinitesimal"` (`ST<:InfinitesimalStrain`), only a
  3-arg `(g,dt,Del)` method exists, so the call throws `MethodError: no method matching
  cpu_oobf_assembly(...)`. The file's own header comment says "only InfinitesimalStrain
  is supported. LogarithmicStrain requires additional handling of the bᵢⱼ left-Cauchy-
  Green tensor" — backwards from what the code actually does. Confirmed pre-existing
  (reproduces identically, modulo `E`- vs `ST`-typed error text, on commit `fbc64ec`,
  before the `AbstractElasticity`-removal refactor that found it) — not introduced by
  that refactor, just newly surfaced by testing a `dynamic_relaxation`+`infinitesimal`
  combination for the first time. `test_collapse.jl` and every other verified
  `dynamic_relaxation` run use the default `strain.deform="finite"`, which is why this
  went unnoticed. Someone needs to either fix `oobf_assembly`'s `ST<:InfinitesimalStrain`
  method to actually get called for infinitesimal runs (the call site would need to stop
  hardcoding `Kc,Gc`), or fix the stale comment and make `"infinitesimal"` genuinely
  unsupported by construction (e.g. an explicit error) rather than an opaque
  `MethodError` deep in a kernel.
- ~~`test/testset/test_collapse.jl` is fully non-functional~~ — **fixed** (branch
  `fix-collapse-problem-pipeline`). `ic_collapse` (renamed `collapse_problem`, matching
  `slump_problem`/`thermal_problem`'s convention) was rewritten to the current
  `get_solver` → `setup_geometry` → `setup_mesh` → `setup_material_constants` →
  `setup_mpts` → `setup_time` → `MechanicalProblem` → `setup_basis` → `export_problem`
  pipeline, the same pattern already used for `collision.jl`/`thermal.jl`. In
  `test_collapse.jl`: `load_simulation_setup` now reads `file["ic"]["problem"]` (a
  `MechanicalProblem`, unpacked into `mesh`/`mpts`/`time`) and `file["cfg"]["solver"]`
  instead of the old `ic["mesh"]`/`ic["mpts"]`/`ic["cmpr"]`/`cfg["instr"]` layout;
  `compute_collapse_error` fixed to call `elastoplasm!(jld2; workflows=[solver])`
  (plural) and replaced the removed `cmpr.ρ0` with `mpts.s.ρ₀[1]` (density is
  per-particle now but uniform across the column); `z0 = copy(setup.mpts.x[end, :])`
  (matrix-style-indexing a `Vector{SVector}`, and reading a `mesh.prprt.dim` field that
  never actually existed on `MeshProperties`) replaced with
  `getindex.(setup.mpts.x, idx)` where `idx = length(setup.mesh.prprt.nel)-1` — the
  same `D`-from-`nel`-length idiom used elsewhere in this codebase (e.g. `slump.jl`'s
  file naming). All 4 `nel` convergence cases (`[5,5]`–`[5,40]`) now pass, including the
  `@test errors[k+1] < errors[k]` monotonic-convergence assertions — a real physics
  correctness signal, not just "doesn't crash." This was also the first real exercise
  of `dynamic_relaxation`/`elastoquasistatic!` this session (the test runs it, not the
  explicit solver) — its previously-**unverified** `cmp`-rerouted reads (see "Typed
  constitutive-model abstraction" below) are now confirmed working.
- **`test/testset/test_workflow.jl`'s full sweep (1 geom × 4 basis × 24
  `strain`/`transfer`/`stab`/`musl` combos = 96 cases) has genuine, reproducible
  numerical-instability failures confined to `smpm` (11/24) and `gimpm` (3/24) —
  `bsmpm`/`mlsmpm` pass 24/24 cleanly.** Earlier notes here claimed these fail
  *silently* via `elastoplasm!(...).success == false` — that's now confirmed **wrong**:
  `elastoplasm`/`elastoplasm!` (`src/home/core/elastoplasm.jl`) have no code path that
  ever sets `success=false`, they return `success=true` unconditionally unless an
  exception propagates out first. Every one of these 14 failures is a real thrown
  exception, caught and logged by the sweep's own `try/@warn` (confirmed by running
  `include("test/testset/test_workflow.jl")` directly instead of through `runtests()`,
  per the guidance below, which surfaces the actual exceptions instead of `runtests()`'s
  swallowed-exception summary). Two distinct exception shapes, and the split lines up
  almost exactly with `transfer.musl`: **12 of 14 failures have `musl=false`** and throw
  `DomainError` (on `locking=true` cases) or `BoundsError: attempt to access N-element
  Vector{SVector{NN,Int64}} at index [...]` (on `locking=false` cases, with the index
  sometimes *negative*) — the `BoundsError` is a particle's topology lookup
  (`p2e`/`p2n`, see `tplgy.jl`) indexing past the mesh entirely, i.e. the particle
  physically left the domain, consistent with an unstable run rather than a data-layout
  bug. The other 2 (both `gimpm`, `locking=true`, `musl=true`) throw the same
  `BoundsError` shape despite MUSL being on, so MUSL reprojection isn't a complete
  explanation — `smpm`+`gimpm` are more failure-prone than `bsmpm`/`mlsmpm` regardless.
  No evidence found of the previously-suspected cross-case state leakage — these look
  like real per-case instabilities, not artifacts of running the sweep in one process.
  **Narrowed (not yet fully root-caused) to the classic MPM cell/grid-crossing
  instability, after an initial wrong theory in this entry was corrected**: the first
  pass here claimed `smpm`'s missing `basis.type` boundary-correction dispatch (which
  `bsmpm` has, see `basis/bsmpm.jl`'s `ϕ∂ϕ`/`t1`-`t4` pieces) was the cause — that's
  wrong. `basis.type` only exists to fix *wide*-stencil kernels (`bsmpm`'s cubic
  B-spline reaches `-1:2`, beyond a particle's own element, so near a domain edge some
  of those nodes don't exist without correction). `smpm`'s stencil is `0:1` — the
  standard bilinear/trilinear hat function using only the particle's own element
  corners — which is compactly supported *within* the mesh by construction and sums to
  1 everywhere, boundary or not, with no correction needed; this matches the existing
  "smpm isolated gradient-consistency glitch" entry below, which found only one
  single-point PoU deviation, not a systemic boundary defect. **The more likely actual
  cause**: `smpm/basis/smpm.jl`'s `N∂N` is a piecewise-linear tent function whose `∂N`
  is discontinuous at `δx=0` (jumps from `1/h` to `-1/h` exactly at a mesh node) — the
  textbook MPM "cell-crossing"/"grid-crossing" instability, where a particle crossing
  an element boundary sees a discontinuous jump in shape-function gradient and thus in
  computed strain-rate/stress at that instant. `bsmpm`'s cubic B-spline is `C¹`-
  continuous by construction (the standard reason B-spline MPM variants exist at all),
  and `mlsmpm`'s per-particle moment-matrix reconstruction likewise avoids the jump —
  both pass 24/24. `gimpm`'s particle-domain averaging (`basis/gimpm.jl`'s `S∂S`)
  smooths the discontinuity partially but not fully, consistent with its much rarer
  3/24 failures vs. `smpm`'s 11/24 with no smoothing at all. This explanation doesn't
  depend on proximity to a domain edge (unlike the boundary-PoU theory), which better
  fits failures spread across the whole `(deform,trsfr)` sweep rather than only
  edge-heavy configurations. Downstream: once a particle's kinematics get jolted by a
  crossing event, `F`/`detF` can drift enough that `volumetric.jl`'s `ΔJp` — which
  raises the F-bar-averaged Jacobian ratio to a fractional power `ΔJ^(1/D)` with no
  sign guard, unlike the zero-mass guard commit `455fc35` added nearby (a real,
  separate un-guarded issue) — receives a negative `ΔJ` and throws `DomainError` under
  `locking=true`; under `locking=false` the drift can instead accumulate until the
  particle's position update pushes it out of the mesh, throwing `BoundsError` in
  `p2e`/`p2n`. `mapsto.jl`'s MUSL branch is basis-kind-agnostic (a generic nodal-
  velocity reprojection/re-solve, unrelated to shape-function continuity), consistent
  with it acting as a general damping/stabilizer that happens to rescue otherwise-
  marginal `smpm` runs rather than interacting with `smpm`'s kernel specifically.
  Most promising next diagnostic step: log `mpts.s.ΔFᵢⱼ`/`detF` (or `basis.∂N`) for a
  single particle across the timestep where it crosses from one element's stencil to
  the next in a failing `smpm` case, and confirm a step-change coincides with the
  crossing. Real fix is likely a higher-order/damped shape function choice for
  affected particles near crossing events (out of scope to design here), or at minimum
  a defensive sign/clamp guard on `ΔJ^dim` in `volumetric.jl` (non-root-cause, but
  closes the crash).
  **Update after `basis.ghost`/`geom.ghost` removal** (mesh padding removed in favour
  of the `iszero(no)` sentinel every kernel already relies on, see the git log for
  that commit): re-ran the `gimpm` 24-case subset in isolation. The composition of
  failures changed but the count didn't meaningfully drop — **5/24 now fail**, all
  `locking=true, musl=false` (`gimpm_finite_tpic`, `gimpm_finite_apic`,
  `gimpm_infinitesimal_std`, `gimpm_infinitesimal_tpic` → `BoundsError`;
  `gimpm_infinitesimal_apic` → `DomainError`). The previously-documented 2
  `locking=true, musl=true` `gimpm` failures **are gone** — consistent with those
  having been driven by the zero-mass ghost nodes the padding introduced (ghost nodes
  never receive particle contributions, so they were a plausible source of the
  `ΔJp`/mass-related instability independent of the grid-crossing theory below). But
  since `musl=false` was already `smpm`'s (and most of `gimpm`'s) dominant failure
  mode before this change, and the grid-crossing/`S∂S`-discontinuity root-cause theory
  above doesn't depend on ghost nodes at all, this confirms ghost-node removal fixed
  one *contributing* failure mode, not the underlying instability — don't read a
  single passing default-parameter `gimpm`+`locking` run (`musl=true` is the default)
  as evidence the bug is fixed; re-run this sweep, not a single case, before revising
  this entry further. `smpm`/`bsmpm`/`mlsmpm` were not re-swept here since their sweep
  config already had `ghost=0` (only `gimpm` used `ghost=1`), so this refactor
  provably can't have changed their behaviour.
  **Re-measured on the `mimic-implementation-logic-of-ample` tip (`61bafbd`), full
  96-case sweep, 2026-08**: the counts quoted throughout this entry are stale. The
  actual current baseline is **smpm 12/24 fail, gimpm 13/24 fail, bsmpm 0/24,
  mlsmpm 0/24 — 25 failures total** (21 `DomainError`, 4 `BoundsError`). The failing
  set is exactly the `musl=false` cases for `smpm` (all 12 of them) plus, for `gimpm`,
  all 8 `locking=true` cases and 2 of the `tpic, lockfalse, muslfalse` ones — the
  gimpm "2 `locktrue, musltrue` failures are gone" claim above is **wrong** at this
  tip; they are back. The typed-strain/stress-tensor port was verified against this
  baseline by running the same sweep on the stashed and unstashed tree in the same
  session: **identical 25-name failure set, identical exception kinds** — i.e. no new
  failures and none accidentally fixed. Full failing case list, for future diffs:
  `smpm_{finite,infinitesimal}_{std,tpic,apic}_lock{true,false}_muslfalse` (12) and
  `gimpm_{finite,infinitesimal}_{std,tpic,apic}_locktrue_musl{true,false}` minus
  `gimpm_finite_tpic_locktrue_musltrue`, plus
  `gimpm_{finite,infinitesimal}_tpic_lockfalse_muslfalse` (13).
- **`ic_collision`'s initial-velocity plot crashes**: `what_plot_field` errors with
  `type NamedTuple has no field data`, because `collision.jl`'s `plot.what` entries are
  hand-built as `(;mpts=(name="u",cblim=(...)))` instead of going through
  `get_mpts_variable_config()[name]` like `slump_problem`'s equivalent plot section does —
  the hand-built `NamedTuple` is missing the `data`/`label`/`unit`/`scale`/`cb` fields
  `what_plot_field` expects. `ic_collision` itself was also just rewritten this session
  (it previously called an undefined `kwargser` function and predated the current
  `get_solver`/`setup_geometry`/`setup_basis` pipeline entirely — that part is fixed),
  but this plot-config bug was pre-existing in the file and wasn't fixed. Workaround:
  run with `plot.status=false` until fixed.

## Planned improvements (not bugs, just known follow-up work)

- **`Problem` type decoupling mesh/point generation from basis — done, as a thin
  wrapper.** `MechanicalProblem{T1,T2,D,E,R,CM,TM,TV,TS} <: AbstractProblem{T1,T2,D,SP}`
  (`src/boot/needs/types/concrete/problem/problem.jl`) bundles `Mesh`+`Point`+`Time` —
  the expensive, IC-defining part of a simulation — as one object, with `SP` tied to
  `PointSolidPhase{T1,T2,D,E,R,CM,TM,TV,TS}` (the type of `mpts.s`), reusing the
  existing `MaterialPointPhase` hierarchy rather than a second, disconnected "phase"
  concept at the problem level, exactly as originally sketched here. `Basis`/`Solver`
  remain separate, un-bundled types (no `Solution{Basis,Solver}`), since `basis.which`
  and `transfer.trsfr` vary independently. The file lives in its own
  `concrete/problem/` subdirectory (not a flat file) specifically so it loads *after*
  `eulerian.jl`/`lagrangian.jl`/`solver.jl` (which define `Mesh`/`Point`/`Time`)
  regardless of alphabetical file order — subdirectories load after all files in their
  parent per the includer's design (see "Project shape"), the same reason
  `concrete/basis/` is a subdirectory rather than a flat file.

  Of the three open questions this entry originally posed: **(1)** `setup_problem`
  (`src/home/init/setup_problem.jl`) is the new entry point — it calls
  `setup_mesh`/`setup_material_constants`/`setup_mpts`/`setup_time` internally and
  returns a `Problem`, replacing that sequence in `slump_problem` (renamed from
  `ic_slump`, since the old name no longer reflected what it returns). `setup_problem`
  no longer takes `mat` as an argument at all — it builds it internally via
  `setup_material_constants(solver)` (renamed from the more cryptic `setup_matconst`; it
  builds the raw scalar elastic/plastic/thermal constants a simulation needs, see its
  docstring in `setup_cmp.jl`), since that function only ever needed `solver` for its
  `T1,T2,D` type parameters, never mesh data — so no caller needs to build a throwaway
  `Mesh`, or build/pass `mat` at all, just to get `setup_problem` going.
  `collapse.jl`/`collision.jl` call `setup_material_constants(solver)` directly the same
  way (they don't use `setup_problem`); `thermal.jl` was left as-is since it's already
  broken independent of this change (calls an undefined `kwargser`, same vintage bug as
  `ic_collapse`/`ic_collision` before their fixes — see Known bugs).
  **(2)** Persistence **was** eventually changed, after initially landing as
  three separate keys: `export_problem` (renamed from `export_setup`, file renamed too,
  matching the `setup_*.jl`/`setup_*` naming convention) now takes
  `(problem::MechanicalProblem, basis, solver, paths; ...)` and writes a single
  `ic["problem"]` entry instead of `ic["mesh"]`/`ic["mpts"]`/`ic["time"]` — see
  "Persistence (JLD2)" below for the full key layout and every reader that had to
  follow (`elastoplasm.jl`, `metanalysis.jl`, `test_performance.jl`, `test_basis.jl`).
  `Basis` stays its own top-level `ic["basis"]` key, unchanged — it's deliberately
  independent of `Problem` (see `MechanicalProblem`'s docstring), so bundling it in
  would misrepresent that relationship. **(3)** Kernel/workflow signatures
  (`workflow!(mpts,mesh,basis,time,solver)`) are **untouched** — `Problem` is, for now,
  only an outer-API/persistence convenience the way option (3)'s first alternative
  described, not threaded through kernels; that remains a larger, currently
  out-of-scope lift if ever wanted. `setup_basis` now takes `(problem::MechanicalProblem,
  solver)` instead of loose `(mesh, mpts, geom, solver)` — the `geom::Geometry` argument
  was dropped entirely (`geom.h` and `mesh.prprt.h` were the same element-size data
  under two names; `setup_basis` now reads `mesh.prprt.h` directly), and `mesh`/`mpts`
  come from `problem` since a `Basis` is always built against a `Problem`'s pair anyway.
  Callers that don't go through `setup_problem` (`collision.jl`) build a
  `MechanicalProblem(mesh, mpts, time)` by hand first — this reordered `time` ahead of
  `basis` there, harmless since `setup_time` doesn't depend on `basis`. The actual
  ergonomic payoff (swap `basis.which` without regenerating `mesh`/`mpts`) works today
  via `setup_basis(problem, solver)`.

  Landing `Problem` also required finishing an abandoned, previously-stashed rename:
  `CartesianMesh`→`AbstractCartesianMesh`, `MeshPhase`→`AbstractMeshPhase`,
  `MaterialPoint`→`AbstractMaterialPoint`, `MaterialPointPhase`→
  `AbstractMaterialPointPhase` (matching the existing `AbstractElasticity`/
  `AbstractRheology`/`AbstractConstitutiveModel` naming convention), and dropping the
  unused `UniformMesh`/`NonUniformMesh` middle layer (`Mesh` now subtypes
  `AbstractCartesianMesh` directly — `NonUniformMesh` had no subtypes anywhere). The
  stash had only renamed `abstract.jl` itself, leaving `eulerian.jl`/`lagrangian.jl`
  inconsistent (still subtyping the old names) — that half was finished as part of this
  work. The stash had also dropped `AbstractSolver`'s `{T1,T2,D}` type parameters,
  which would have broken ~20 call sites dispatching on `AbstractSolver{T1,T2,D}` —
  that part of the stash was **not** kept; `AbstractSolver{T1<:Integer,T2<:Real,D}`
  stays parameterized.
- **Kernel granularity**: `elast!` is the model to follow — decomposed into small
  functions each doing one specific task (e.g. computing log strain), rather than one
  monolithic kernel body. Other large kernels (e.g. `update.jl`'s `elastoplast`/
  `elasto` dispatch functions, which branch on `strain.deform`/`basis.how` inline)
  would benefit from being split the same way.
- **Typed constitutive-model abstraction — done.** `mpts.s.cmp::Vector{CM} where
  CM<:AbstractConstitutiveModel` (`src/boot/needs/types/concrete/constitutive.jl`) now
  bundles the *static* elastic+plastic material constants (`Gc`, `Kc`, `Del`, `Hp`,
  `c₀`, `cᵣ`, `ϕ₀`) into one typed object per particle, replacing the old global
  `cmpr::NamedTuple` that was threaded through ~15 function signatures and persisted as
  `ic["cmpr"]` — that key and the whole `cmpr` argument are gone. Evolving plastic state
  (`Δλ`, `ϵpII`, `ϵpV`) deliberately stays *outside* `cmp`, as flat mutable vectors
  directly on `PointSolidPhase` (unchanged), since it's written every timestep and
  bundling it into an immutable per-particle struct would force a full-struct
  reconstruction on every write. This mirrors, but does not copy, the
  `AbstractConstitutiveModel`/`PerfectlyElastic`/`DruckerPrager` idea prototyped on
  branch `fancy-implementation` — re-implemented against the current architecture
  rather than merging that branch's incompatible `Point`/`Basis`/`D`-as-type-param
  shape backward, and resolving that branch's `Gc`/`Kc` scalar-vs-vector inconsistency
  by always using one convention: one `CM` instance per particle, plain scalar `T2`
  fields inside it (see `constitutive.jl`'s docstrings for the full rationale,
  including why `Del` is sized `nstr×nstr`, the Voigt stiffness matrix, not `D×D` as
  originally sketched). Two things worth knowing if you touch this next:
  - `PerfectlyElastic` is defined but has **no construction path** — `setup_cmp`
    (`src/home/init/setup_cmp.jl`) always builds `DruckerPrager`, since both live
    retmap kernels (`DP.jl`, `J2.jl`) need its full field set. It's currently dead
    scaffolding, the same category `PointSolidPhase.elast::E` used to be in before it
    was removed (see below).
  - The old duplicate flat `PointSolidPhase.c₀`/`cᵣ`/`ϕ₀` vectors have been retired —
    `DP.jl`/`J2.jl` already read `cmp[p].c₀`/`.cᵣ`/`.ϕ₀` instead, and the plotting path
    (`get_coh0`/`get_phi0` in `src/home/plot/variables.jl`) was the last reader,
    rerouted to `[cmp.c₀ for cmp in mpts.s.cmp]`/`[cmp.ϕ₀ for cmp in mpts.s.cmp]`. The
    fields are gone from the struct and its constructor in `setup_mpts.jl`.
  - `PointSolidPhase.rheo::R`/`DruckerPragerRheology`/`AbstractRheology` — the
    parallel, never-read component-field attempt at per-particle constitutive data
    (built every `setup_mpts` call, alongside `elast`, but consumed nowhere) — have
    been retired entirely, since `cmp` already fully supersedes what `rheo` was
    reaching for. Removed: the `AbstractRheology{T1,T2}` abstract type
    (`abstract.jl`), the `DruckerPragerRheology` struct (`lagrangian.jl`), the `R`
    type parameter from `PointSolidPhase`/`Point`/`MechanicalProblem`, and the `rheo`
    field/its construction in `setup_mpts.jl` — plus the now-stale `R`/
    `R<:AbstractRheology` from every kernel signature that pattern-matched on
    `Point{T1,T2,...,E,R,...}` (~20 files across `basis/`, `core/common/`,
    `core/solver/explicit/`, `core/solver/dynamic_relaxation/`, `plot/`). At the time
    `rheo` was removed, `elast::E` was deliberately left alone — it was separately-
    flagged dead scaffolding (see next bullet), not the same "duplicate of `cmp`"
    problem `rheo` was. It has since been removed too, once the tensor port gave `E`
    an actual redundant replacement (`ST`) — see next bullet.
  - **`tensor.jl`'s typed `AbstractStrain`/`AbstractStress` abstraction — done**, see
    "Typed strain/stress tensor storage" immediately below.
  - **`dynamic_relaxation`'s rerouted reads — verified for the plain path, still
    unverified for u-P.** `implicit.jl` (`pt_solve!`, `pt_solve_uP!`, `relax`,
    `mapsto_uP`, `elastoquasistaticuP!`, `elastouP!`) was mechanically converted from
    `cmpr.Kc/Gc/c` to `mpts.s.cmp[p]` the same way the `explicit` path was.
    `elastoquasistatic!`/`relax`/`pt_solve!` (the plain, non-u-P path) are now confirmed
    working — exercised for real by `test_collapse.jl`'s convergence sweep (branch
    `fix-collapse-problem-pipeline`), including passing its `@test errors[k+1] <
    errors[k]` correctness assertions. The u-P variant (`pt_solve_uP!`/`mapsto_uP`/
    `elastoquasistaticuP!`/`elastouP!`) remains unverified — and per the separate
    known-bug entry above, throws immediately if exercised regardless (`fint_p2n!`
    dispatch never registered), so its `cmp` rerouting specifically hasn't been proven
    correct either way.
- **Typed strain/stress tensor storage — done.** `mpts.s.ϵᵢⱼ`/`.ϵn`/`.σᵢ`/`.σn`/`.τᵢ` no
  longer hold bare `SMatrix{D,D}`/`SVector{nstr}` values. They hold typed tensor objects
  from `src/boot/needs/types/concrete/tensor.jl`, each storing a volumetric+deviatoric
  additive split:
  - `LogarithmicStrain{S,T,L}` / `InfinitesimalStrain{S,T,L}` (`vol::T`,
    `dev::SMatrix{S,S,T,L}`) — `ϵᵢⱼ`/`ϵn`, picked by `solver.strain.deform`
    (`"finite"`→log, `"infinitesimal"`→small strain) in `build_solid_phase`.
  - `CauchyStress{S,T,L}` / `KirchhoffStress{S,T,L}` (`p::T` positive in compression,
    `dev::SMatrix{S,S,T,L}`) — `σᵢ`/`σn` are always Cauchy, `τᵢ` always Kirchhoff.

  The abstract supertypes `AbstractTensor{T}`/`AbstractStrain{S,T,L}`/
  `AbstractStress{S,T,L}` live in `types/abstract.jl`, **not** in `concrete/tensor.jl` —
  `lagrangian.jl` constrains `PointSolidPhase`'s type parameters against them and loads
  first under the alphabetical includer, so declaring them in `tensor.jl` would fail to
  load. `PointSolidPhase`/`Point`/`MechanicalProblem` gained three **trailing** type
  parameters `ST<:AbstractStrain, SC<:AbstractStress, SK<:AbstractStress` after `TS`;
  trailing specifically so the ~40 kernel signatures matching only `Point{T1,T2,D}` /
  `Point{T1,T2,D,E}` needed no edit. Only `setup_mpts.jl` and `problem.jl` list the full
  parameter set. `TS` (the Voigt `SVector` type) is now unused by
  `PointSolidPhase`'s fields but was kept as a parameter to avoid churn.

  Ported from `fancy-implementation`, but merged into the current architecture rather
  than merging that branch's `Point`/`Basis` shape backward. What is a straight port:
  `LogarithmicStrain`/`KirchhoffStress`, `get_tensor`/`get_voigt`, `Base.getindex` on
  strains, the `LinearAlgebra.eigen` overload, `_trial_elastic_strain`,
  `_trial_elastic_stress`. What is **new code** (that branch left
  `InfinitesimalStrain`/`CauchyStress` as bare structs with no methods): both types'
  constructors, `_infinitesimal_strain`, the full-Voigt-vector constructors on both
  stress types, `Base.zero` for all four, and the `get_J2`/`get_τII` invariants.

  Four things worth knowing before touching this:
  - **The stored split is representational, not canonical.** The only invariant the
    types guarantee is that `get_voigt`/`get_tensor` reconstruct the tensor. `dev` is
    *not* required to be trace-free in `S` dimensions, because two different pressure
    conventions genuinely coexist in this codebase and both predate the port:
    `_trial_elastic_stress` uses `tr(ϵ)/3` even in 2D (plane strain — the out-of-plane
    component is real physics, just never stored), while the Voigt constructors use
    `-tr(σ)/S`, matching `σTr` (`retmap/DP.jl`) and `yield_J2` (`retmap/J2.jl`), which
    divide by 2 in 2D. `get_J2`/`get_τII` therefore re-derive a trace-free deviator from
    `get_voigt` rather than trusting the stored `dev`, so they agree whichever writer
    produced the object. **Always consume a stored tensor via `get_voigt`/`get_tensor`,
    never by reading `.dev`/`.p` and assuming they mean what you'd expect.**
  - **`get_J2` was NOT unified with `retmap/J2.jl`'s same-named helper** (which the
    original plan called for). They compute different things: the retmap one returns
    `(‖ξ‖, n̂)` — the deviator norm and yield-surface normal — while the tensor one
    returns the scalar invariant `J₂`. Merging them would either change what a call site
    receives or force the return mapping onto the non-canonical stored split. The retmap
    helper was **renamed `_yield_normal`** instead, so the collision is gone.
  - **Return-mapping kernels stay in Voigt space internally.** `DP.jl`/`J2.jl` read
    `get_voigt(...)` once, run the unchanged closed-form/CPA algebra on `SVector`s, and
    wrap the result once at the store. Same for the Jaumann-rate `elast`/`elast_fast`
    kernels. This is deliberate: it keeps the numerics bit-identical (see below) and
    preserves `elast_fast`'s hand-inlined performance, while the *stored* value is still
    genuinely typed. `drucker_prager` goes one better — `σn` already produces a pressure
    and a scaled deviator separately, so the returned `KirchhoffStress` is built from
    those directly instead of being re-split from a Voigt vector.
  - **Numerical impact, measured, not assumed.** The Kirchhoff-stress predictor and the
    Voigt algebra are bit-identical (20000/20000 random cases, max rel diff exactly
    0.0). What is *not* bit-identical is the storage round-trip itself: splitting a
    tensor into `(p,dev)` and reassembling costs ≤1 ulp of the largest component
    (measured: exactly 1.0 ulp worst case over 20000 random 2D and 3D stresses;
    ~54%/51% of cases are still bit-exact). Through an `eigen`-based log-strain update
    that ≤1-ulp input perturbation grows to ≤6.8e-13 (2D) / 5.0e-14 (3D) relative, and
    through the DP return map to ≤1.8e-13 relative on `Δλ`/`ϵpII`. `test_workflow.jl`'s
    96-case sweep produced a **byte-identical failure set** before and after (25
    failures, same 25 case names — see Known bugs), and `test_performance.jl` reports
    **-12.1% memory / -5.4% allocs / 0.0% time** overall vs the stored baseline. Every
    new type's construction, `get_voigt`, `get_tensor`, `get_J2`, `get_τII`,
    `_trial_elastic_stress`, `_infinitesimal_strain` and `zero(...)` measure **0 bytes**
    allocated in a 10 000-iteration loop inside a concrete function. (Beware: measuring
    these one-shot through a generic `probe(name, f, args...)` helper reports a spurious
    16–112 B floor from the harness's own dynamic dispatch — the pre-port `SVector`
    baselines show the identical floor. Use loop-accumulating probes.)
  - JLD2 needs **no** special handling: the structs are plain immutables over `T` and
    `SMatrix`. Verified by an explicit save/load round-trip of a full `Point` in both
    `deform` modes — types and values identical, `eltype`s preserved.
- ~~`PointSolidPhase.elast::E` (component field typed `AbstractElasticity`) remains an
  intentionally unused dead field~~ — **removed**, along with `AbstractElasticity`,
  `FiniteElasticity`, and `LinearElasticity` themselves, and the `E` type parameter from
  `Point`/`PointSolidPhase`/`MechanicalProblem`. This entry originally argued `elast`
  should be left alone pending evidence it duplicated something — that evidence arrived
  once the typed strain/stress tensor port landed (see "Typed strain/stress tensor
  storage" above): `E<:FiniteElasticity`/`E<:LinearElasticity` was, by that point, a
  pure dispatch tag distinguishing finite vs. infinitesimal strain, and `ST` (`
  LogarithmicStrain`/`InfinitesimalStrain`) already encodes exactly that distinction —
  confirmed by grep that `mpts.s.elast`/`.elast` had zero reads anywhere in `src/`/
  `test/` before removal, same zero-hit result the original entry already reported.
  Every kernel that used to dispatch on `E<:FiniteElasticity`/`E<:LinearElasticity` —
  `elast.jl`'s `elast`/`elast_fast` (2 methods each) and `dynamic_relaxation/
  {fint,update}.jl` (2 methods each) — now dispatches on `ST<:LogarithmicStrain`/
  `ST<:InfinitesimalStrain` instead, via the full positional parameter list up to `ST`
  (`Point{T1,T2,D,CM,TM,TV,TS,ST}`), since Julia can't skip earlier type parameters to
  bind a later one. ~10 more files that only used `Point{T1,T2,D,E}` as an
  unconstrained 4-parameter placeholder (not real dispatch) had that free variable
  renamed, no logic change. Verified: clean load (no silent include failures), full
  `plast.constitutive ∈ {"DP","J2"} × strain.deform ∈ {"finite","infinitesimal"}` smoke
  matrix plus `nonloc.status=true` and `perf.status=true` (`elast_fast` path) all
  `success=true`, `dynamic_relaxation`/`elastoquasistatic!` passes with
  `strain.deform="finite"` (the only combination it ever supported — `"infinitesimal"`
  throws `MethodError` in `fint.jl`'s `oobf_assembly`, confirmed pre-existing and
  unrelated to this change by reproducing the identical error on the prior commit),
  `test_workflow.jl`'s 96-case sweep unchanged at 25/96 failures (same per-basis-kind
  distribution: smpm 12, gimpm 13, bsmpm 0, mlsmpm 0), `test_collapse.jl` 4/4,
  `test_basis.jl` passing, and `test_performance.jl` **-17.2% memory / -5.4% allocs /
  0.0% time** vs the stored baseline — better than the tensor port's own -12.1%, as
  expected since this also stops building `elast`'s dead zero-filled arrays every
  `setup_mpts` call.
- ~~`PointSolidPhase.P::Vector{T2}` field~~ — **removed**, same pass as `elast::E`
  above. Confirmed zero reads anywhere in `src/`/`test/` before removal; pressure lives
  on the stress objects themselves now (`mpts.s.σᵢ[p].p`/`mpts.s.τᵢ[p].p`, or
  `get_voigt`'d), making the separate scalar field redundant rather than merely
  unused.
- **`get_vector` renamed `get_voigt`, and its strain-side shear convention was wrong —
  fixed, and `mutate`/`_mutate` retired as a result.** The user asked whether the old
  free functions `mutate(ϵ,Χ,type)`/`_mutate(ϵ)` (`elast.jl`, predating the tensor port)
  duplicated what `tensor.jl` already provided. They didn't — `get_vector(strain::
  AbstractStrain)` used the **tensor**-shear convention (`γxy=εxy`, no doubling), but
  every place strain interacts with the elastic stiffness matrix `Del` needs the
  classical **engineering**-shear Voigt convention (`γxy=2εxy`), which is *why* `Del` is
  built the way it is (`σ=Del·ε_voigt` only works with the ordinary isotropic stiffness
  matrix when `ε_voigt` carries that factor of 2) — `mutate`'s `Χ=2.0`/`Χ=0.5` branches
  existed purely to patch that mismatch at every `Del`-facing call site. The user's key
  insight: the *engineering*-shear vector is the actual "conventional" Voigt strain
  vector — `get_vector`'s "no factor 2" convention was the anomaly, not `mutate`'s
  patch. Confirmed **zero blast radius** before fixing it: `get_vector` was never once
  called on a strain object anywhere in `src/` (grep-confirmed), so correcting the
  convention changed no live numerics, only fixed an as-yet-unexercised part of the API
  before anything came to depend on the wrong version. Renamed `get_vector`→`get_voigt`
  throughout (one name, dispatched by multiple dispatch onto `AbstractStrain`/
  `AbstractStress` methods — not two names; only the `AbstractStrain` methods' behavior
  changed, `AbstractStress` has no engineering/tensor ambiguity and was renamed only).
  Added `LogarithmicStrain(ϵ::SVector{3/6,T})`/`InfinitesimalStrain(ϵ::SVector{3/6,T})`
  full-engineering-Voigt-vector constructors (`tensor.jl`), mirroring the pattern
  `CauchyStress(σ::SVector)`/`KirchhoffStress(τ::SVector)` already used on the stress
  side — `vol=(ϵ[1]+ϵ[2])/3` (2D, matching the existing plane-strain-3D convention, NOT
  stress's `tr/S`), `dev` built by halving shear into tensor form. Exact inverse of
  `get_voigt`: `get_voigt(LogarithmicStrain(ϵ)) == ϵ` (verified to machine epsilon,
  ~7e-18, over 20000 random 2D/3D cases).

  With that gap closed, 8 of `mutate`'s 10 call sites collapse into existing or new
  `tensor.jl` calls: the two `Χ=1.0,:tensor` stress-reshape sites (`elast.jl`,
  `dynamic_relaxation/fint.jl`) become `get_tensor(mpts.s.σᵢ[p])` directly; the two
  `Χ=2.0,:voigt` strain→`Del` sites become `get_voigt(InfinitesimalStrain(ϵ))`; the two
  `Χ=0.5,:tensor` `Del\\...`→strain sites (`DP.jl`, `J2.jl`) become
  `LogarithmicStrain(cmp.Del\\τᵢ)` directly, dropping the intermediate `SMatrix` wrap
  entirely. The two remaining sites (`elast.jl`/`fint.jl`'s Jaumann-rate correction term
  `σJ*ω'+σJ'*ω` — genuinely not "the" particle's stress, so no typed object exists for
  it) got a small dedicated replacement, `voigt_of(M)` (`elast.jl`, stress-convention
  only, no `Χ` parameter needed since only one case remains). `_mutate` (already
  zero-caller dead code) and `mutate` were both deleted outright.
  `retmap/MCRetMap.jl`'s own `mutate(...,0.5,:tensor)` call (line 44) was **left
  untouched** — that file is dead/unwired (only referenced in a commented-out line,
  `update.jl:47`, a known separately-documented gap), so it now calls an undefined
  function if ever exercised; not in scope to fix alongside this.

  Verified: clean load; round-trip identity machine-epsilon exact; every replaced call
  site numerically equivalent to the old `mutate`-based code to machine epsilon or
  bit-exact (0.0 for the stress-side and Jaumann-scratch sites); zero allocations on
  every new/changed function; full `DP`/`J2` × `finite`/`infinitesimal` smoke matrix
  plus `nonloc`/`perf_fast` all `success=true`; the `fint.jl` infinitesimal kernel
  (unreachable through `elastoquasistatic!` per the known bug below) verified directly
  by calling `oobf_assembly!` on a hand-built `Point`/`Mesh`/`Basis` — ran to completion
  with finite, non-NaN output; `dynamic_relaxation` with `strain.deform="finite"` still
  passes; `test_workflow.jl`'s 96-case sweep unchanged (71 passed/25 failed, same total
  as both prior baselines — this refactor touches only basis-kind-agnostic per-particle
  bookkeeping, so a total match here is conclusive); `test_collapse.jl` 4/4,
  `test_basis.jl` passing; `test_performance.jl` unchanged (-17.2% memory/-5.4%
  allocs/0.0% time vs baseline — this refactor's touched code isn't on the benchmarked
  finite-strain default path, so no shift was expected).
- **3D conformity check**: verify 3D simulations behave consistently across the
  explicit solver's configuration space (basis kind, `stab.locking`, `strain.deform`)
  — most of this session's verification was 2D-only.
- **`basis.how`/`basis.ghost` are GIMP-specific concepts living in the generic `basis`
  config section.** `how` (particle-domain update mode) is read unconditionally in
  `update.jl`'s dispatch regardless of which basis kind is actually active, and
  `ghost` similarly sits at the top level even though it only meaningfully matters for
  `GimpBasis`'s wider stencil. Worth exploring moving both to live on `GimpBasis`
  itself (construction-time fields, defaults, and validation co-located with the kind
  that actually uses them) rather than as basis-kind-agnostic top-level knobs that
  `bsmpm`/`smpm`/`mlsmpm` silently ignore.
- **`Basis.N`/`∂N` storage rework**: currently plain `Matrix{T2}`/`Array{T2,3}`
  specifically to dodge the StaticArrays allocation-elision limits documented below —
  correct and fast, but a workaround rather than a considered data-layout design; worth
  revisiting once/if a cleaner approach is found that doesn't reintroduce the
  allocation regression.
- **`smpm`'s isolated gradient-consistency glitch**: `test_basis.jl` found
  `max|Σᵢ∂Nᵢ|` off by exactly `1.0` at one sweep point for `LinearBasis`, unlike
  `gimpm`'s broadly-degraded boundary PoU. Never root-caused — an error of exactly
  `1.0` smells like a real (if narrow, single-point) bug rather than float noise, worth
  a closer look.
- **`DynamicRelaxationSolver`'s fate**: now that `solution` picks between
  `ExplicitSolver`/`ImplicitSolver` only (see "Controlling solver behaviour" below),
  this third struct is even more clearly orphaned than before. Decide: wire it in as a
  third `solution` value with real semantics, or state plainly (here, or in
  `solver.jl`) that it's intentionally dormant scaffolding, so nobody "completes" the
  selection logic to include it without knowing whether that's actually wanted.
- **Staged path to a poro-hydro-mechanical solution: mechanical solver (already
  exists) → standalone fluid-only (hydro) solver → fuse the two.** Intent is to build
  and validate the hydro solver independently before attempting coupling, rather than
  designing the coupled solve first. The fluid-only solver needs `PointFluidPhase`
  (`lagrangian.jl`) to gain real fields — it's currently a literally empty placeholder
  struct (`# Add concrete fields as needed`), `mpts.f` is always `nothing` in
  practice — plus its own workflow entry point, e.g. a `hydrodynamic!` in
  `explicit/worflow.jl` mirroring `elastodynamic!`, with its own P2G/G2P/update kernel
  path (mirroring `mapsto`/`update`'s structure) before any
  solid-fluid coupling work is meaningful.
- **Thermal solution — fixed and now viable, on branch `fix-thermal-phase-shape-function-cache`.**
  `thermal_problem` (renamed from `ic_thermal`, `src/home/script/example/thermal.jl`,
  rewritten to mirror `slump_problem`'s `get_solver → setup_geometry → setup_problem →
  setup_basis → export_problem` pipeline) → `elastoplasm(jld2;
  workflows=[thermodynamic!])` now runs end-to-end and produces physically sane output
  (verified: temperature stayed bounded within `[boundary, initial]` and cooled
  monotonically from the boundary inward, both with `musl=true` and `musl=false`).
  What was actually broken, beyond the previously-known `mpts.ϕ∂ϕ` reads:
  - `mpts.ϕ∂ϕ[nn,p,:]` reads in `std_p2n`'s `MeshThermalPhase` overloads
    (`mapsto/fun/p2n/std.jl`), `augm_p2n` (`mapsto/fun/augm.jl`), `picflip_n2p`
    (`mapsto/fun/n2p.jl`), and `heatflux` (`update/fun/heatflux.jl`) all replaced with
    `basis.N[nn,p]`/`∂Nrow(basis.∂N,nn,p,Val(D))`, the same pattern the solid phase
    already used.
  - `std_p2n`'s 3D thermal overload was mistyped as `mesh::Mesh{T1,T2,3}` instead of
    `mesh::MeshThermalPhase{T1,T2,3}` — never actually dispatched as the thermal
    kernel it clearly was (body reads `mpts.t`/`mesh.cᵢ`/`mesh.mcT`/`mesh.oobq`).
  - `heatflux`'s signature never bound `D` (needed for `Val(D)`); `thermo`'s
    (`update/update.jl`) signature bound a stale `E` positionally into `Point`'s
    now-`D` slot post-`rheo`-removal. Both fixed.
  - **The real missing piece**: `thermodynamic!` (`explicit/worflow.jl`) calls
    `mapsto(mpts, mesh.t, basis, dt, instr)`, but only a solid-phase `mapsto` existed
    (different arg count/order, `mesh::Mesh` not `MeshThermalPhase`, requires a
    gravity vector). A new `mapsto(mpts::Point,mesh::MeshThermalPhase,basis,dt,solver)`
    overload was added in `mapsto.jl`, mirroring the solid one's p2n→solve→n2p(→musl
    augm reprojection) structure but without gravity/APIC/finite-strain concerns —
    every individual kernel it calls (thermal `std_p2n`/`euler`/`picflip_n2p`/
    `augm_p2n`/`augm_solve`) already existed and needed no new design, just this one
    orchestrating wrapper.
  - `mesh.t`/`mpts.t` were unconditionally `nothing` — `setup_mesh`/`setup_mpts`/
    `setup_problem` all gained a `thermal::Bool=false` kwarg (default `false`, so
    every other caller is unaffected) that builds a real zero-initialized
    `MeshThermalPhase`/`PointThermalPhase` instead. Only `thermal_problem` sets
    `thermal=true`.
  - `get_thermal.jl`'s initial temperature was `zeros(nmp).*mat[:initial_temperature]`
    (always zero regardless of config) — fixed to `ones(nmp).*...`.
  - **Not fixed, left as a known gap**: `thermal_problem` defaults to
    `plot.status=false`. `ic_thermal`'s old plotting config used a stale
    string-tuple format (`what=[("mpts","T"),("mesh","T")]`) incompatible with the
    current `get_plot_field`/`get_mpts_variable_config()`-based API, and even the
    modern API has no `"T"` entry in `get_mpts_variable_config()`/
    `get_mesh_variable_config()` yet, plus `display.jl` carries two *different*,
    incompatible `what_plot_field(mesh::Mesh,...)` dispatch styles (an old
    string-`opts.what`-based one that reads `mesh.t.T` directly, and the new
    NamedTuple-config style `get_plot_field` actually calls) — unifying these is a
    separate, real piece of follow-up work, out of scope for getting the solve loop
    itself working.
  - **Separately discovered, fixed in the same branch**: 3D APIC transfer for the
    *solid* phase (`mapsto/fun/p2n/apic.jl`) had the identical dead-`ϕ∂ϕ` bug, plus a
    mismatched `mesh::MeshSolidPhase` type (every other `p2n!` kernel is called with
    the top-level `mesh::Mesh`) and inconsistent `mesh.m`/`mesh.mv` vs `mesh.s.m`/
    `mesh.s.mv` field access. Fixed to match the working 2D APIC kernel's pattern.
    Verified the fix itself is correct (700+ iterations of real physics before
    failure), but 3D+APIC then hits a separate, pre-existing `DomainError` (sqrt of a
    negative value) — the same general explicit-MPM numerical-instability class
    already documented elsewhere in this file for other undertested basis/config
    combinations — persisting even with `stab.locking=false`. Root cause not chased
    further; 3D+APIC was already flagged as untested territory under "3D conformity
    check" below, and 1D APIC has the same `mesh.m`/`mesh.s.m` inconsistency,
    unfixed (out of scope — 1D is essentially never exercised in this codebase).
- **Add a JLD2 time-series export option, as an alternative to plot-only checkpoints.**
  Every workflow's time loop (`explicit/worflow.jl`) computes `checks =
  sort(collect(time.t[1]:solver.plot.freq:time.te))` and calls `bake(mpts,mesh,t,solver)`
  at each checkpoint — but `bake` (`src/home/plot/bake.jl`) only does anything if
  `solver.plot.status` is true, and in that case it unconditionally renders/saves a PNG
  via `get_plot_field`. There is currently no way to record a field's evolution as raw
  data: to get numbers out today you must either turn on plotting and get pixels, not
  values, or write ad hoc postprocessing against final-state-only JLD2 output. A natural
  fix reuses the same `checks` cadence and `solver.plot.what`-style field selection to
  instead (or additionally) append selected `mpts`/`mesh` field values into a JLD2
  dataset at each checkpoint — effectively a `plot.status`-style toggle (e.g.
  `export.status`/`export.what`) sitting next to `bake()` in the workflow loop rather
  than replacing it. Writing the export is only half the feature: nothing today reads
  a saved time series back out either — `metanalysis.jl`'s `extract_field_data`/
  `postprocess_fields` is the closest existing precedent (it already loads a field via
  `get_mpts_variable_config()` and plots it, but only ever from *final-state* JLD2
  output across many single-timestep simulation runs, never a time axis within one
  run). A companion loader/plotter that opens one run's exported time-series dataset
  and plots a field's evolution over `t` (reusing `get_mpts_variable_config()` for
  field metadata the same way `metanalysis.jl`/`bake.jl` already do) would need to be
  built alongside the export itself, not as a separate follow-up — an export nobody
  can read back is not actually useful on its own.

## Precision (32-bit) and StaticArrays performance gotchas

`dtype=(;T0=(Int32,Float32),bits=Int32(32),...)` is a supported but rarely-exercised
path — the default is always `(Int64,Float64)`, so bugs that only manifest under
`T1=Int32`/`T2=Float32` silently pass every day-to-day test. Two concrete patterns
found and fixed by actually running an `Int32`/`Float32` `slump_problem`→`elastoplasm`
end-to-end (not just skimming the code):

- **`@index(Global)` inside a `@kernel` function is always `Int64`**, regardless of the
  solver's index type `T1`. Passing it straight into anything dispatching on `::T1`
  (e.g. `eval_basis(mpts,mesh,basis,ip::T1,nn::T1)`) throws a `MethodError` under
  `T1=Int32` — silently fine under the default `T1=Int64` only by coincidence. Always
  wrap it: `T1(p)`.
- **Bare numeric literals (`2.0`, `0.5`, `1/3`, ...) silently promote to `Float64`**,
  even inside a function generic over `T2`. Innocuous under the default `T2=Float64`,
  but under `T2=Float32` it produces a mix of `Float32`/`Float64` values that eventually
  hits a `MethodError` at whatever call site requires them to match (often several
  functions downstream of where the literal actually is — the error site is not the bug
  site). Always write `T2(2.0)` etc., never a bare literal, in any function generic over
  a float type parameter.

Separately, when caching per-particle results into a `Vector{SVector}`/`Vector{SMatrix}`
(or building one from scratch) inside a hot loop: **bundling `NN` computed values into
one whole-object `SVector{NN,...}`/`SMatrix{NN,D,...}` can silently heap-allocate**,
even when every intermediate step is individually zero-alloc in isolation — Julia's
escape analysis has a complexity/size cutoff past which it gives up eliding the
allocation, and this was empirically true regardless of *how* the object was built
(`MVector` buffer, `ntuple`/`Val`, `hcat`, a flat tuple — all ~allocated the same
amount for `NN=16`). Two mitigations that measurably worked here (see `Basis.N`/`∂N`
above): (1) `ntuple(f, n)` needs `n` as `Val(n)`, not a plain `Int`, to compile-time
unroll — even when `n` is itself a type parameter, if it's used as a runtime value it
won't unroll; (2) store the *whole cached collection* as a plain `Matrix`/`Array` (like
`Point`'s pre-existing `Δnp`/`Dᵢⱼ` fields) written via scalar element assignment,
rather than as a `Vector` of per-particle `SVector`/`SMatrix` objects — reading a row
back out via `M[nn,:]` still allocates for a runtime `nn`, so read it back via a small
`Val`-unrolled `ntuple` (see `∂Nrow` in `basis.jl`) instead of array-slicing syntax.
Always verify with `@allocated`/`@code_typed` in a wrapping *function* (not top-level
REPL scope, which has its own unrelated allocation noise) rather than assuming a
StaticArrays-based kernel is allocation-free because it "looks" that way.

## Persistence (JLD2)

Simulation setup is saved/loaded as `ic["problem"]` (a `MechanicalProblem`, bundling
`mesh`/`mpts`/`time`) and `ic["basis"]` (kept as its own top-level key — deliberately
independent of `Problem`, see `MechanicalProblem`'s docstring); solver config as
`cfg["solver"]` (NOT `cfg["instr"]` — that key name is stale/wrong despite appearing
in some scripts). `ic["mesh"]`/`ic["mpts"]`/`ic["time"]` as three separate keys is an
**old layout, no longer written** — `export_problem` writes `ic["problem"]` only now.
`elastoplasm`/`elastoplasm!` unpack `problem = file["ic/problem"];
mesh,mpts,time = problem.mesh, problem.mpts, problem.time; basis = file["ic/basis"]`;
`elastoplasm!`'s write-back after each workflow rebuilds
`file["ic/problem"] = MechanicalProblem(mesh,mpts,time)` (mesh/mpts are mutated
in-place by the workflow, so this just rewraps the same objects) alongside
`file["ic/basis"] = basis`. Both `elastoplasm`/`elastoplasm!` take a
`workflows::Vector{Function}` kwarg (plural — not `workflow`). **`ic["cmpr"]` no
longer exists** — constitutive material constants now live per-particle on
`mpts.s.cmp` (persisted as part of `ic["problem"].mpts`), see "Typed constitutive-model
abstraction" under Planned improvements. Other readers of the old three-key layout were
updated to match: `metanalysis.jl` (`reference["ic/problem"].mesh`/`.mpts`),
`test/testset/test_performance.jl` (`ic["problem"].mesh`/`.mpts`/`.time`), and
`test/testset/test_basis.jl` (`load(jld2, "ic/problem")`). `test_collapse.jl` was
**not** updated — it was already fully non-functional before this change (see Known
bugs: undefined `kwargser`, stale `cfg["instr"]`, matrix-style `mpts.x` indexing, a
nonexistent `ic["cmpr"]` key already predating this rename) and remains so; fixing it
is one more thing to account for whenever someone actually picks that file up.

## Operating the package

Typical end-to-end run (this is also the standard smoke test after any refactor):

```julia
using ElastoPlasm
L,nel = [64.1584, 64.1584/4.0], [40,10]      # domain size, element counts (2D here; 3 entries = 3D)
jld2  = slump_problem(L, nel; cli()...)           # builds mesh/mpts/basis/time, saves setup, returns the .jld2 path
out   = elastoplasm(jld2; workflows = [elastodynamic!, elastoplastic!])
@assert out.success
```

- `slump_problem(L, nel; fid="...", kwargs...)` (`src/home/script/example/slump.jl`) is the
  reference example problem. It calls, in order: `get_solver` → `setup_geometry` →
  `setup_problem` (internally: `setup_mesh` → `setup_material_constants` (takes
  `solver` directly, not a `Mesh` — see `setup_problem`'s docstring and the `Problem`
  planned-improvement entry) → `setup_mpts`, which itself calls `setup_cmp` to build the
  per-particle `mpts.s.cmp` → `setup_time`; returns a `MechanicalProblem`) →
  `setup_basis`, then `export_problem(...)` to persist everything to a `.jld2` file and
  returns that path.
  Use it as the template for defining a new example problem.
- `cli()` (`src/home/api/solver/cli.jl`) parses interactive/CLI overrides for solver
  config; `cli(; ui=true)` prompts interactively, `cli()` with no args picks defaults
  non-interactively — pass its result as `kwargs...` into `slump_problem`/`get_solver`.
- `get_solver(; dim=2, kwargs...)` (`src/home/api/solver/get_solver.jl`) merges
  user kwargs over `get_default()`, resolves the execution backend, and builds the
  `Cairn` dispatch table (`cairn.ignite`/`.mapsto`/`.update`/`.implicit`) — this is
  where solver-type kernels actually get instantiated (`SomeKernel(CPU())`).
- `elastoplasm(jld2_path; workflows=[...])` reopens the `.jld2`, unpacks
  `(mesh, mpts, basis, time, solver)`, and runs each workflow function in order
  as `workflow!(mpts, mesh, basis, time, solver)`. `elastoplasm!` is the in-place
  variant (opens the file `"r+"` so postprocessing can write back into it). Both accept
  any callable matching that signature as a "workflow" — the built-ins are
  `elastodynamic!`/`elastoplastic!` (explicit solver, `explicit/worflow.jl`) and
  `elastoquasistatic!` (dynamic_relaxation solver, `dynamic_relaxation/worflow.jl`);
  they can be freely composed/chained in the `workflows` vector since they all share
  the same `(mpts, mesh, basis, time, solver)` call signature.
- Plotting is opt-in via `solver.plot.status` (set through `slump_problem`'s `plot=(...)`
  kwarg or `cli()`); when on, `elastoplasm` auto-saves a PNG per workflow named
  `"$(dim)_$(basis.which)_$(solvertype)_$(deform)_$(workflow)_$(quantity).png"` under
  the run's `dump/.../plot` path (`basis.which` was added so runs that only differ by
  basis kind don't silently overwrite each other's output).

### Running tests / benchmarks

```julia
julia --project=. -e 'using Pkg; Pkg.test()'          # or: include("test/runtests.jl")
```

`test/runtests.jl` interactively lets you pick which `test/testset/test_*.jl` files to
run (auto-runs all of them under `GITHUB_ACTIONS=true`). Note it swallows exceptions
inside `runtests()`'s `try/catch` (just increments a fail counter) — when debugging a
failing test file, `include()` it directly instead of through `runtests()` to see the
real stacktrace, e.g.:

```julia
using Test,JLD2,ProgressMeter,Suppressor,Plots,LaTeXStrings,REPL.TerminalMenus,ElastoPlasm
using BenchmarkTools, KernelAbstractions
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync
global ROOT, DATASET = joinpath(pwd(),"test"), joinpath(pwd(),"test","dataset")
include("test/testset/test_performance.jl")
```

`test_performance.jl` benchmarks core kernels directly via `solver.cairn.*` (including
`ignite.shpfun!`) and compares against `test/dataset/performance_baseline.jld2` (keyed
by CPU name); delete/update that file deliberately if a change is a genuine, intended
perf shift rather than a regression.

`test_basis.jl` builds a small `n×m` mesh and, for every basis kind, sweeps a particle
along the mid-height node row, plotting each tracked node's `Nᵢ(x)`/`∂Nᵢ/∂x(x)` and the
partition-of-unity `Σᵢ Nᵢ(x)` — useful both as a quick correctness check on a new/changed
basis kind and as a way to visually confirm findings like the `GimpBasis` boundary PoU
issue (see Known bugs). Both this file and `test_performance.jl` report their numeric
checks informationally (`println`, not hard `@test` gates) — deliberately, since not
every basis kind is expected to hold e.g. exact PoU everywhere, so a strict `@test`
would be a false failure rather than a real regression signal. When adding LaTeX to a
generated figure (`LaTeXStrings.jl`, already a project dependency): use it for genuine
math expressions (axis labels, legend entries like `L"N_i(x)"`) but keep plain
identifiers — basis kind names, etc. — as plain strings; wrapping a multi-letter plain
word in math mode (e.g. `\mathrm{bsmpm}`) does not render reliably under the default GR
backend's mathtext support and was an actual bug hit while building `test_basis.jl`.

### Controlling solver behaviour via `defaults.jl`

`get_default()` (`src/home/api/solver/defaults.jl`) returns just `default`, the base
`NamedTuple` config merged by `get_solver` — this is the single place that defines
every tunable knob and its out-of-the-box value. (It used to also return a hardcoded
`solver_type` — `ExplicitSolver` always, even for `dynamic_relaxation` runs — that's
gone now; see `solution` below for how the concrete type is chosen instead.)

- `solution` — `"explicit"`/`"implicit"`, plain user-facing string picking the concrete
  solver struct at construction time in `get_solver.jl` (`ExplicitSolver`/
  `ImplicitSolver` — both share the same field layout, `DynamicRelaxationSolver` is
  unused scaffolding, not wired into this selection). This is also what
  `elastoplasm.jl`'s generated filenames and `logs.jl`'s startup banner print
  (`solver.solution`) — previously they used `nameof(typeof(solver))`, which printed
  the literal type name `"ExplicitSolver"` even for `dynamic_relaxation` runs, since
  every run constructed that one type regardless of what it was actually doing.
  Picking `solution="explicit"` and then calling an implicit-only workflow (or vice
  versa) now fails with a `MethodError` at the workflow call, not silently — those
  workflow entry points (`elastodynamic!`/`elastoplastic!` in `worflow.jl`) are typed
  to their specific concrete solver struct on purpose. The shared `ignite()`
  (`core/common/ignite.jl`) is the one place dispatching on the abstract
  `AbstractSolver` instead, since both `elastodynamic!`/`elastoplastic!` *and*
  `elastoquasistatic!` (`dynamic_relaxation`) call it regardless of `solution`.
- `dtype` — arithmetic precision (`bits`, element types `T0`)
- `basis` — `which` (`"bsmpm"`/`"gimpm"`/`"smpm"`/`"mlsmpm"`, see `get_basis`), `how`
  (GIMP domain update mode: `"detFij"`/`"Fii"`/`"Uii"` finite, `"detΔFij"`/`"ΔFii"`/
  `"ΔUii"` infinitesimal), `ghost` (ghost-node count)
- `strain` — `deform` (`"finite"`/`"infinitesimal"`)
- `transfer` — `trsfr` (transfer scheme: `"std"`/`"tpic"`/`"apic"`), `C_pf` (PIC/FLIP
  blend), `musl` (MUSL velocity reprojection on/off)
- `stab` — `locking` (F-bar volumetric locking correction on/off), `damping`

  (`strain`/`transfer`/`stab` used to be one flat `fwrk` block — split by concern since
  `deform`, the transfer-scheme knobs, and the stabilization knobs are three genuinely
  separate things that happened to share a vague name. `ExplicitSolver`/`ImplicitSolver`
  now carry three corresponding fields instead of one `fwrk` field.)
- `bcs` — `dirichlet` boundary condition matrix, one `[lower upper]` row per dimension
- `grf` — Gaussian random field generator for heterogeneous cohesion/friction fields
  (`status` toggles it on; see `GRF.jl`)
- `plast` — `status`, `constitutive` (`"DP"`/`"J2"`/`"MC"`/`"camC"` — not all are
  wired up, check `update/update.jl`'s `init_update` dispatch before relying on one)
- `nonloc` — non-local plastic strain regularization (`status`, `ls` length scale)
- `plot` — `status`, `freq` (plot every N `Time` checkpoints), `dpi`, `what` (list of
  field specs to plot, keyed by name via `get_mpts_variable_config()`)
- `perf` — `status`; when true, forces `deform="infinitesimal"` and disables `nonloc`
  and swaps in the `_fast` kernel variants (see `init_update`) for a lighter-weight run
- `backend` — `select` (execution backend, `"host"`/GPU target) and `distributed`

Two ways to change solver behaviour:

1. **Per-run override (preferred, no code change)** — pass any of these keys as
   keyword arguments through `slump_problem`/`get_solver`; `get_solver` only keeps kwargs
   whose keys already exist in `default` (unrecognized kwargs just get a `@warn`, not
   an error), then deep-merges them over `default` via `merge(ref, user)`. Note this is
   a shallow-per-key merge — to override one field of a nested `NamedTuple` block
   (e.g. just `transfer.trsfr`) you must pass the *whole* `transfer = (; trsfr=...,
   C_pf=..., musl=...)` NamedTuple, not just the one field, since `merge` replaces the
   whole `:transfer` entry rather than recursing into it. Example:
   ```julia
   slump_problem(L, nel; basis=(;which="gimpm",how="Uii",ghost=0), strain=(;deform="finite"), transfer=(;trsfr="apic",C_pf=1.0,musl=true), stab=(;locking=true,damping=0.1))
   ```
2. **Change the package-wide default** — edit the literal value in `get_default()`
   directly. Do this only for a genuine change of the shipped default behaviour, not
   for one-off experimentation (use kwargs for that) — every example/script that
   doesn't override a given key inherits from here.

New solver-config knobs must be added as a new key inside the relevant block (or a new
top-level block) in `get_default()`'s `default` NamedTuple, then threaded through
wherever `init_ignite`/`init_mapsto`/`init_update`/`init_implicit` (in `get_solver.jl`
and the `explicit`/`dynamic_relaxation` solver files) branch on `instr[:section][:key]`
to select kernels — grep those `init_*` functions for the existing pattern
(`instr[:transfer][:trsfr] == "apic"`-style branches) before adding a new one. If the
new field belongs on the solver struct itself (like `solution` did), it also needs
threading through `ExplicitSolver`'s/`ImplicitSolver`'s field list
(`src/boot/needs/types/concrete/solver.jl`) and the corresponding positional arg in
`get_solver.jl`'s final constructor call — both structs are kept in lockstep by hand,
there's no shared base struct to edit once.

### Using `cli()`

`cli(; ui::Bool=false)` (`src/home/api/solver/cli.jl`) turns `get_default()`'s config
tree into a `Dict{Any,Any}` of kwargs suitable for splatting into `slump_problem`/
`get_solver` as `cli()...`. Two modes:

- `cli()` (default, `ui=false`) — **non-interactive**, just returns `get_default()`'s
  values reshaped into the kwargs `Dict`. This is what you want for scripts and
  automated runs (`slump_problem(L, nel; cli()...)`) — it's the plain-defaults path, not a
  prompt.
- `cli(ui=true)` — **interactive**: first shows a `MultiSelectMenu` letting you pick
  which top-level sections to configure (`solution`, `dtype`, `basis`, `strain`,
  `transfer`, `stab`, `grf`, `plast`, `nonloc`, `plot`, `perf`, `backend`); for each
  selected section it walks
  `get_option()` (the parallel tree of `("prompt", [choices])` tuples) and shows a
  `RadioMenu`/`MultiSelectMenu` per leaf option (nested sections like `grf.param.Iₓ`
  are handled recursively; sections with a `status` field skip their sub-prompts when
  you answer `false`; `plot.what` gets a `MultiSelectMenu` instead of `RadioMenu` since
  multiple fields can be plotted at once). Unselected sections silently fall back to
  their `get_default()` value, not a prompt. Use this only from a real terminal/REPL —
  it blocks on `TerminalMenus.request`, so never call `cli(ui=true)` from a
  non-interactive script or from within an agent session.
- To override specific options while still using `cli()` for everything else, splat
  `cli()` first then override afterward — later kwargs win in Julia's `merge`-by-
  keyword-splat semantics: `slump_problem(L, nel; cli()..., basis=(;which="gimpm",
  how="Uii", ghost=0))`. Same shallow-merge caveat as `get_default()` applies: you must
  pass the full nested `NamedTuple` for a section, not just the one field you want to
  change.
- `get_option()` is also useful standalone as a reference for what values are actually
  valid per key (e.g. `get_option().transfer.trsfr` → `["std", "tpic", "apic"]`) — treat it
  as the source of truth for valid choices, since some of `get_default()`'s tunables
  (e.g. `plast.constitutive`) accept values `get_option()` doesn't fully enumerate
  (check `init_update`'s dispatch in `update/update.jl` if in doubt whether a value is
  actually wired up).

## Conventions / house rules

- `NN` (nodes-per-element count) belongs to `Basis`, never to `Point` or `Mesh` — it's
  topology, not point/mesh state.
- Prefer `SVector`/`SMatrix` (StaticArrays.jl) over runtime-sized arrays for per-point
  fixed-size data.
- When converting one solver path to a new convention, apply the identical mechanical
  pattern already validated on the sibling path (explicit ↔ dynamic_relaxation) rather
  than reinventing it — these two solver directories are meant to stay structurally
  parallel.
- Verify refactors by actually running the `slump_problem` → `elastoplasm(...; workflows=[...])`
  pipeline end-to-end, not just a clean-load check — a file can `include` cleanly and
  still be functionally broken (wrong kwarg name, stale dict key, matrix-style indexing
  into what's now a `Vector{SVector}`, etc.).
- `test/testset/test_performance.jl` benchmarks core kernels directly via
  `solver.cairn.*`; keep it in sync with the `Cairn` dispatch-table shape whenever that
  shape changes.
- `using ElastoPlasm` flushes (deletes) everything under `dump/` on load
  (`rootflush`, `src/boot/needs/utils.jl`) if that directory already has contents.
  Never run more than one `julia --project=. -e '...'` invocation against this repo
  concurrently while any of them is still writing to `dump/` — a second session's
  `using ElastoPlasm` will delete the first session's still-in-progress output out from
  under it (this produced a very confusing "No such file or directory" `savefig` error
  mid-session before the actual cause — concurrent sessions, not a real bug — was found).

## Commit message conventions

This repo follows Conventional Commits: `<type>: <Summary>`, summary capitalized, no
trailing period. Pick the type by what actually changed, not by what the change was
*for*:

- `feat:` — new functionality or capability that wasn't there before.
- `fix:` — corrects behavior that was wrong (a bug, a crash, an incorrect result).
- `refactor:` — restructures existing code without changing external behavior (e.g.
  renaming, moving logic between files, simplifying a dispatch pattern).
- `chore:` — maintenance with no production/behavior impact: build tooling, dependency
  bumps, config, repo housekeeping, formatting.
- `docs:` — documentation-only changes (README, CLAUDE.md, docstrings).
- `test:` — adding or updating tests only, no source changes.
- `perf:` — a performance improvement, not a correctness fix.

`git log --oneline` in this repo also has some older/looser types in its history
(`misc:`, `bug:`) — treat those as historical, not part of the current convention;
prefer the list above for new commits.
