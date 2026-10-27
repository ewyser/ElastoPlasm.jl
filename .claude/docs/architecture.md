## Core types

- `Mesh{T1,T2,D}` — Eulerian background grid. Carries no connectivity (`e2n`/`e2e` live
  on `Basis`).
- `Point{T1,T2,D,CM<:AbstractConstitutiveModel,ST<:AbstractStrain,EL<:AbstractElasticLaw,L}`
  — material points (Lagrangian). Every parameter but `L` is a real dispatch axis; `L == D*D`
  only exists to keep the `D×D` kinematic fields (`SMatrix{D,D,T2,L}`: `F`, `ΔF`, `∇v`, ...)
  concrete, since Julia can't compute `D*D` in a field type —
  kernels never dispatch on it and recover a field's type via `eltype`. `ST` is the typed
  strain storage on `mpts.s` (see "Typed strain/stress tensor storage" in `planned-improvements.md`);
  `mpts.s.σᵢ[p]` returns a `CauchyStress`, **not** an `SVector` — read it with
  `get_tensor(...)`, or its `.p`/`.dev` split. Stress and strain are full 3D
  tensors in every dimension (plane strain in 2D), stored as pressure/volumetric plus a trace-free
  Voigt deviator, while kinematic fields (`F`, `ΔF`, `∇v`, `ω`) stay `D×D`; see
  "Typed strain/stress tensor storage" in `planned-improvements.md`. Carries no `NN` or
  connectivity (those live on `Basis`).
  `mpts.x :: Vector{SVector{D,T2}}` — NOT a matrix; index with `getindex.(x, i)` or
  iterate, never `x[i, :]`. `mpts.s.m` is the constant solid mass `(1−n₀)ρ₀Ω₀`, set once
  at setup; kernels read it instead of recomputing `ρ·Ω`, which `deform!` keeps equal to it
  (a future fluid phase needs its own, separately evolving mass). `mpts.s.cmp::Vector{CM}` is the per-particle
  constitutive-model bundle (moduli and plastic constants only; the elastic laws use `Kc`/`Gc`
  directly; `get_Del(cmp)` in `elasticity.jl` builds the Voigt stiffness when a matrix is needed) (see "Typed constitutive-model abstraction" in
  `planned-improvements.md`). `CM` resolves to `DruckerPrager` or `VonMises` depending
  on `material.plastic` (see "DP/J2 retmap kernel unification"). `ST` (`LogarithmicStrain`/
  `InfinitesimalStrain`) and `EL<:AbstractElasticLaw` (`HenckySolid`/`ImprovedHenckySolid`/
  `HypoelasticSolid`, the elastic-law tag `elast.jl` and `get_dt`'s `_Ktan` dispatch on) are both
  picked together from `material.elastic` by `build_solid_phase`.
- `Problem{T1,T2,D,P,B,S}` (`src/boot/needs/types/problem/problem.jl`) — what is solved: `mpts`,
  `mesh`, `basis`, `time`, `solver`, built by `Problem(case)`; what workflows run on
  (`elastodynamic!(problem)`). The `Case` it was built from is not a field (see "`Problem`" in
  `planned-improvements.md`).
- `Simulation{P}` (`src/boot/needs/types/problem/simulation.jl`) — a problem ready to run:
  `problem`, `workflow` (a name, one of `WORKFLOWS`), `lstps`, `folder` (absolute; a name given is
  under `dump/`), built by `Simulation(problem; workflow, lstps, folder)`; the default workflow follows the
  solver (`"elastoquasistatic!"` for an `ImplicitSolver`), and an explicit workflow on an implicit solver is
  refused there rather than as a `MethodError` at `execute!`. A `.jld2` is a saved one
  (`save_simulation`/`load_simulation`); `execute!(sim)` / `execute!(jld2)` run it in place (`execute` on a copy). The chain is
  `Case` → `Problem` → `Simulation`, each built from the one above.
- `Case` (`src/boot/needs/types/problem/case.jl`) — the plain-data description a problem is built from: domain, zones
  (`Zone` = shape `Polygon`/`Block`/`Ball`/`Region` + material overrides), default material
  (`MATERIAL`), boundary conditions, solver overrides, thermal flag, run settings. `case_points`
  turns it into `PointData`, the explicit per-point contract `setup_mpts` reads. One build path
  (`setup_simulation`), one TOML format (`save_case`/`load_case`); the drawing window edits it. See
  "Describing a case" in `operations.md`.
- `Basis{T1,T2,D,NN,K<:AbstractBasis,TR<:AbstractTransfer}`
  (`src/boot/needs/types/basis/basis.jl`) — independent struct owning ALL
  mesh/point connectivity (`e2n`, `e2e`, `p2n`, `p2e` — `e2e` also doubles as the
  neighbor-search structure for nonlocal plastic-strain regularization, see
  `.claude/bug/fixed/nonlocal-regularization-on2-and-asymmetry-bug.md`), plus `kind::K`
  (dispatches shape-function evaluation, `eval_basis`/`shpfun!`), `transfer::TR`
  (dispatches P2G/G2P transfer-scheme kernels, `p2n!`/`Bij` — see "Transfer scheme
  dispatch" below), `NN` (nodes per element — only known once mesh+basis kind are fixed,
  hence its own type param), `type` (per-axis node boundary-layer classification,
  `BSplineBasis`-specific), and the shape-function cache `N`/`∂N` (see below).
- Construction order is always **mesh → material points → basis** (basis depends on
  both). Function argument order is always **`(mpts, mesh, basis, ...)`**.
- Shape-function evaluation is `eval_basis(mpts, mesh, basis, ip, nn)`. Concrete kinds
  live in `basis/{bsmpm,gimpm,smpm,mlsmpm}.jl`, each defining its own
  `stencils::NTuple{D,UnitRange{T}}` field on the kind struct.
- **Shape values/gradients are cached, not recomputed per call site.** `Basis.N` /
  `Basis.∂N` (plain `Matrix{T2}` / `Array{T2,3}`, sized `(NN,nmp)` / `(NN,D,nmp)` — see
  `gotchas.md` for why not `Vector{SVector}`/`Vector{SMatrix}`)
  are populated once per particle per timestep by the `shpfun!` ignite kernel
  (`src/home/core/common/shpfun.jl`). Every P2G/G2P/update kernel reads `basis.N[nn,p]` /
  `∂Nrow(basis.∂N,nn,p,Val(D))` (a small helper in `basis.jl` extracting a gradient row
  via a `Val`-unrolled `ntuple`, since `M[nn,:]` on an `SMatrix` allocates for a runtime
  `nn`) instead of calling `eval_basis` inline. `dynamic_relaxation` does the same: its
  shape functions are cached by `ignite` once per load step, which is valid because
  particles only move after the step has converged.
- `MLSBasis` (`basis/mlsmpm.jl`) is a fourth kind implementing structured-grid Moving
  Least Squares shape functions (Cao et al. 2025, *Comput. Mech.* 75:655–678, §2.3.1 —
  only the structured-grid case). Needs no boundary node-type correction (`basis.type`)
  — a per-particle moment-matrix projection restores partition-of-unity everywhere on
  its own — but needs a once-per-particle `(D+1)×(D+1)` matrix build+inversion, which is
  why the `shpfun!` caching above exists. `basis.which` and `basis.trsfr` are fully
  orthogonal — any basis kind pairs with any transfer scheme, including `mlsmpm`+`std`.

**Transfer scheme dispatch via `Basis`'s `TR` type parameter.** `p2n!` (P2G projection)
and `Bij` (APIC's affine-velocity update) are dispatched on `Basis`'s `TR`
(`StdTransfer`/`TpicTransfer`/`ApicTransfer`), the exact same "config string → concrete
marker instance → multiple dispatch" pattern `kind`/`get_basis` already uses for basis
kind, and the same pattern the DP/J2 `retmap` unification uses for `Point`'s `CM`+`ST`
(see "DP/J2 retmap kernel unification" in `planned-improvements.md`).
`get_transfer(basis.trsfr, T2, D, nmp)` (`transfer.jl`) maps `"std"`/`"tpic"`/`"apic"`
to `StdTransfer()`/`TpicTransfer()`/`ApicTransfer{T2,D,L}(nmp)`; `init_mapsto` registers
`p2n!(CPU())`/`Bij(CPU())` unconditionally, no string branch. `Bij`'s
`TR<:StdTransfer`/`TpicTransfer` methods are no-ops, so `mapsto.jl` calls it
unconditionally too. Unrecognized `basis.trsfr` strings fail fast in `get_transfer` at
`setup_basis` time.

`ApicTransfer{T,D,L}` carries `Bᵢⱼ`/`Dᵢⱼ` (`Vector{SMatrix{D,D,T,L}}` — `L` is a separate
type parameter, computed as `D*D` only once `D` is concrete inside the inner
constructor, since `D*D` on a bare `TypeVar` isn't valid in a struct's own field-type
declaration). These fields used to live on `Point` unconditionally (`Δnp`/`Bᵢⱼ`/`Dᵢⱼ`,
100% APIC-exclusive, allocated on every particle even under `std`/`tpic`); `Δnp` was
dead code and got deleted, `Bᵢⱼ`/`Dᵢⱼ` moved onto `ApicTransfer` (`basis.transfer.Bᵢⱼ`,
not `mpts.Bᵢⱼ`). They persist with the basis (`problem.basis`), not on `mpts`. Removing
this dead/wasted per-`Point` allocation measured **-18.3% memory / -5.4% allocs / 0.0%
time** on `test_performance.jl` vs. the pre-refactor baseline.

**Config note**: `basis.trsfr`/`.C_pf` (transfer scheme + PIC/FLIP blend) live under the
`basis` config section (mirroring `Basis` owning both `kind` and `transfer`); `stab.musl`
(MUSL velocity reprojection) lives under `stab`, since it's a stabilization technique
applied regardless of which transfer scheme ran, not transfer-scheme-specific. There is
no `transfer` config section anymore. See `operations.md` for the full config layout and
the per-field merge (a shallow merge caused a real bug twice during this reorganization).

**Known, deliberate behavior change**: thermal-only workflows now work under any
`basis.trsfr` (previously threw `MethodError` for `apic`/`tpic`) — that restriction was
an accident of `std_p2n`/`tpic_p2n`/`apic_p2n` being separately-named kernels (thermal
only ever had a method under the `std_p2n` name); now that every `p2n!` method shares
one name, thermal's method (which carries no `TR` constraint, since it never needed one)
matches regardless of `TR`. Confirmed via `git worktree` diff against the pre-refactor
commit to be the *only* behavior change from this refactor; not a regression.

## Explicit vs dynamic_relaxation solvers

- `src/home/core/common/` — kernels/wrappers shared by **both** solver paths:
  `tplgy.jl` (`p2e2n`, topology), `shpfun.jl` (`shpfun!`/`Dij_nd`), `ignite.jl`
  (`init_ignite`/`ignite()`), `domain.jl` (`particle_domain(basis.kind, F, ℓ₀, ℓ)`, GIMP's particle-domain
  update, a no-op for the other kinds: called per particle where `F` is updated, in the explicit `deform` and
  the dynamic relaxation's end-of-load-step `update`; no kernel or `basis.which` branch of its own). Which
  update is the `GimpBasis`'s `domain` field (`StretchDomain`/`UndeformedDomain`, from `basis.domain` by
  `get_domain`, `types/basis/gimpm.jl`), dispatched on in `particle_domain`; another update is one
  `AbstractDomain` type and one method. Both `elastodynamic!`/`elastoplastic!` *and*
  `elastoquasistatic!` call the same `ignite(mpts,mesh,basis,solver::AbstractSolver)`
  wrapper, dispatching on the *abstract* solver type since `dynamic_relaxation` doesn't
  require a particular concrete `solution` value (unlike `elastodynamic!`/
  `elastoplastic!`, typed to `ExplicitSolver` specifically — see `solution` in
  `operations.md`).
- `src/home/core/solver/explicit/` — the primary, actively-used path
  (`elastodynamic!`, `elastoplastic!`). Kernel dispatch table lives on
  `solver.cairn` (a `NamedTuple`, dot-accessed: `solver.cairn.ignite.tplgy!`,
  `solver.cairn.mapsto.map.p2n!`, `solver.cairn.update.deform!`, ...), built by
  `init_ignite`/`init_mapsto`/`init_update`/`init_implicit` in `get_solver.jl`. Every
  kernel/function in this path takes `basis::Basis` as its 3rd positional arg after
  `(mpts, mesh, ...)`.
- `src/home/core/solver/dynamic_relaxation/` — implicit quasi-static path
  (`elastoquasistatic!`), the dynamic-relaxation (DR) MPM of Sutcliffe et al. (2026,
  IJNME 127:e70384), Algorithm 1 without aggregation. Gravity is applied in `lstps`
  uniform load steps (keyword, default 10, up to 9.81); each step is solved by `relax!`
  (`implicit.jl`) with nodal `Δu` as the primary variable, held in `mesh.s.v`. One DR
  iteration is `dr_kinematics` (ΔF, F, and `ϵᵢ ← ϵn` so no state carries between
  iterations) → the explicit path's `elast!` → `dr_fint` (Eqs 19-22); then the Underwood
  critical damping (Eq 40, capped at 1.9) and the midpoint update (Eqs 32-33, with Eq 44's
  `2+αΔT` denominator; Eq 32's `2+ΔT` is a typo). The fictitious mass is the diagonal
  Hencky P-wave estimate of Eqs 35/38/50 (`fictitious_mass`, `ζ_M=2`, no `N` weight),
  assembled once per step where `ΔF=I`. Converges on `‖R‖/‖fᵉˣᵗ‖ < tol` (1e-9) over free
  DOFs; the workflow returns the iteration count per load step. `mpts.s.u` accumulates the total displacement over the workflow (the explicit path stores the per-step increment `dt·v` there instead). Finite strain only
  (`hypoelastic` throws). Pass options through `execute!` as a closure:
  `workflows=[(a...)->elastoquasistatic!(a...; lstps=40)]`.

  **Eq 35 takes one element size `h`**, so on elongated cells (`h = minimum(h)`) the
  long direction is over-massed and the Rayleigh-quotient damping stalls near
  `‖R‖/‖fᵉˣᵗ‖≈1e-8` for tens of thousands of iterations. Use square cells, as the paper
  does. On the paper's §4.1 column (gimpm, one element wide) it reproduces the Fig 5
  trend (stress error 1.3e-2 → 3.8e-4 for h = 6.25 → 0.78 with 40 load steps) at
  ~100–1500 iterations per load step.

  `plastic=true` runs `retmap!` after `elast!` in every iteration, restoring the step's
  converged `ϵpII` (with `ϵpII[2] ← ϵpII[1]`, as the explicit path seeds it) first;
  `G` sets the final gravity. On the §4.3 von Mises collapse (2² MPs/cell, 10 load steps) the
  run-out is 5.74 m against the paper's ~5.6–5.8 m (Fig 17a), **reading the paper's `ρ_y` as a
  Lode radius** (`√(2J₂) = ρ_y`, i.e. this repo's `κ = c₀ = ρ_y`), not the `√J₂ = ρ_y` Eq 49
  prints: with `√J₂` the block barely slumps (2.2 m), since a smooth-base heap can stand up to
  `H = 4c_u/γ`, which is then exactly the 8 m block height. Convergence is far slower than
  the paper's (1.2M iterations; one load step hits the 1e6 cap), and the cause is identified as
  **small cuts**: the slow mode lives on nodes just beyond the free face whose lumped mass is
  1e-5–1e-4 of the maximum but whose fictitious mass is 12–25% of it (Eq 35 has no `N`
  weight), so they are over-massed ~10⁴× relative to their stiffness. The paper's own fix is
  aggregation; this repo fixes nodes with a low nodal solid fraction instead (`stab.sfrac`,
  `common/sfrac.jl`, DR only — see "Nodal solid fraction" in `planned-improvements.md`), which
  brings this run to ~16k iterations at a <1% run-out bias. Re-assembling `M̃` every iteration (§4.3's variants, with
  Eq 51's `M_p = K` for yielding points) was tried and did not help, so the mass stays
  assembled once per load step.

## Open/resolved issues

Currently-open bugs live under `.claude/bug/known/`, one file per issue. Resolved bugs with
their full investigation history live under `.claude/bug/fixed/`. Check both before
re-investigating a symptom from scratch.
