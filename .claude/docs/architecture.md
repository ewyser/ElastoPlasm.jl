## Core types

- `Mesh{T1,T2,D}` — Eulerian background grid. Carries no connectivity (`e2n`/`e2e` live
  on `Basis`).
- `Point{T1,T2,D,CM<:AbstractConstitutiveModel,ST<:AbstractStrain,EL<:AbstractElasticLaw,L}`
  — material points (Lagrangian). Every parameter but `L` is a real dispatch axis; `L == D*D`
  only exists to keep the `D×D` tensor fields (`SMatrix{D,D,T2,L}`, `CauchyStress{D,T2,L}`,
  `KirchhoffStress{D,T2,L}`) concrete, since Julia can't compute `D*D` in a field type —
  kernels never dispatch on it and recover a field's type via `eltype`. `ST` is the typed
  strain storage on `mpts.s` (see "Typed strain/stress tensor storage" in `planned-improvements.md`);
  `mpts.s.σᵢⱼ[p]` returns a `CauchyStress`, **not** an `SVector` — read it with
  `get_voigt(...)`. Carries no `NN` or connectivity (those live on `Basis`).
  `mpts.x :: Vector{SVector{D,T2}}` — NOT a matrix; index with `getindex.(x, i)` or
  iterate, never `x[i, :]`. `mpts.s.m` is the constant solid mass `(1−n₀)ρ₀Ω₀`, set once
  at setup; kernels read it instead of recomputing `ρ·Ω`, which `deform!` keeps equal to it
  (a future fluid phase needs its own, separately evolving mass). `mpts.s.cmp::Vector{CM}` is the per-particle
  constitutive-model bundle (see "Typed constitutive-model abstraction" in
  `planned-improvements.md`). `CM` resolves to `DruckerPrager` or `VonMises` depending
  on `material.plastic` (see "DP/J2 retmap kernel unification"). `ST` (`LogarithmicStrain`/
  `InfinitesimalStrain`) and `EL<:AbstractElasticLaw` (`HenckySolid`/`ImprovedHenckySolid`/
  `HypoelasticSolid`, the elastic-law tag `elast.jl` and `get_dt`'s `_Ktan` dispatch on) are both
  picked together from `material.elastic` by `build_solid_phase`.
- `MechanicalProblem{T1,T2,D,CM,ST,EL,L} <: AbstractProblem{T1,T2,D,SP}`
  (`src/boot/needs/types/problem/problem.jl`) — bundles `mesh::Mesh`+
  `mpts::Point`+`time::Time` as the IC-defining part of a simulation, built via
  `setup_problem` (see "`Problem` type decoupling..." in `planned-improvements.md`).
  Deliberately excludes `Basis`/`Solver` — those stay independent types.
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
  `nn`) instead of calling `eval_basis` inline. `dynamic_relaxation`'s ~14 call sites are
  the one exception — no ignite phase there, so they call `eval_basis` directly.
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
not `mpts.Bᵢⱼ`). They persist as part of `ic["basis"]` now, not `ic["problem"]`. Removing
this dead/wasted per-`Point` allocation measured **-18.3% memory / -5.4% allocs / 0.0%
time** on `test_performance.jl` vs. the pre-refactor baseline.

**Config note**: `basis.trsfr`/`.C_pf` (transfer scheme + PIC/FLIP blend) live under the
`basis` config section (mirroring `Basis` owning both `kind` and `transfer`); `stab.musl`
(MUSL velocity reprojection) lives under `stab`, since it's a stabilization technique
applied regardless of which transfer scheme ran, not transfer-scheme-specific. There is
no `transfer` config section anymore. See `operations.md` for the full config layout and
a real shallow-merge bug this reorganization surfaced twice.

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
  (`init_ignite`/`ignite()`). Both `elastodynamic!`/`elastoplastic!` *and*
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
- `src/home/core/solver/dynamic_relaxation/` — quasi-static path (`elastoquasistatic!`,
  plus the not-fully-wired u-P variant `elastoquasistaticuP!`/`elastouP!`). Runs
  end-to-end through the same `slump_problem`/`elastoplasm(jld2; workflows=[...])`
  pipeline as the explicit solver — but is much slower to smoke-test (a full
  `elastoquasistatic!` run takes ~1 minute vs. seconds for `elastodynamic!`/
  `elastoplastic!`), so budget for that when verifying a change touches this path.

## Open/resolved issues

Currently-open bugs live under `.claude/bug/known/`, one file per issue. Resolved bugs with
their full investigation history live under `.claude/bug/fixed/`. Check both before
re-investigating a symptom from scratch.
