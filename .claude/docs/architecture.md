## Core types

What is solved (a `Case`, EasyDraw's: a layout with its settings) and how (a `Solver`) meet in the `Problem`, which a
`Simulation` runs. EasyDraw describes a model for any solver and gives its values no meaning; ElastoPlasm reads them
(`case_physics`, `Solver(case)`). The solver holds nothing of a case, so one solver solves several cases of its
dimension; the case's physics reaches the problem only through `Problem(case, solver)`, and after setup the problem's
own types carry it (`laws(mpts)`).

```
EasyDraw: what a model is, for any solver; plain data, what a file holds
┌─ Layout: where things are ───────────────────────────────┐   ┌─ Settings: the values, by zone name ─────────────────────┐
│ Layout(; L, nel, ni, zones)                              │   │ Settings(; zones, defaults, sides, options, run)         │
│                                                          │   │                                                          │
│ L      ::Vector{Float64}  size [m]                       │   │ zones    ::Dict{String,Dict{String,Any}}                 │
│ nel    ::Vector{Int}      cells per axis                 │   │ defaults ::Dict{String,Any}                              │
│ ni     ::Int              points per cell, axis          │   │ sides    ::Dict{String,String}                           │
│ zones  ::Vector{Zone}     a name, a shape each           │   │ options  ::Dict{String,Any}  physics, numerics           │
│ folder ::String           its file's folder              │   │ run      ::Dict{String,Any}  stages, workflow, folder    │
└──────────────────────────────────────────────────────────┘   └──────────────────────────────────────────────────────────┘
                             │                                                              │
                             └──────────────────────────────┬───────────────────────────────┘
                                                            ▼ Case(layout, settings)
                               ┌─ Case: a layout with its settings ───────────────────────┐
                               │ Case(layout, settings=Settings())  copies, checked       │
                               │                                                          │
                               │ layout   ::Layout                                        │
                               │ settings ::Settings                                      │
                               │                                                          │
                               │ save_case, load_case; zone_at, points, partition, …      │
                               │ draw(case; catalog), assign(case; catalog): the window   │
                               └──────────────────────────────────────────────────────────┘

ElastoPlasm: how a case is read and solved
┌─ Reading a case: functions, no type ─────────────────────┐   ┌─ Solver: how it is solved ───────────────────────────────┐
│ case_settings(; properties, material, bcs, solver,       │   │ Solver(case; run, overrides…)  reads case.settings       │
│                 thermal, run) -> Settings                │   │ Solver(; dim, sections…)                                 │
│ case_physics(case) -> (; material, properties, bcs,      │   │                                                          │
│                          laws, grf, thermal)             │   │ solution ::String      explicit | implicit               │
│ case_folder(case), case_points(case, mesh)               │   │ dtype    ::NamedTuple  precision: T1, T2                 │
└──────────────────────────────────────────────────────────┘   │ basis    ::NamedTuple  which, trsfr, C_pf, domain        │
                                                               │ stab     ::NamedTuple  locking, damping, musl, sfrac     │
                                                               │ nonloc   ::NamedTuple  status, ls                        │
                                                               │ plot     ::NamedTuple                                    │
                                                               │ output   ::NamedTuple                                    │
                                                               │ backend  ::NamedTuple  select, exec                      │
                                                               │ run      ::NamedTuple  te, tg, tep, workflow, lstps      │
                                                               │ cairn    ::NamedTuple  the kernels                       │
                                                               └──────────────────────────────────────────────────────────┘
                             │                                                              │
                             └──────────────────────────────┬───────────────────────────────┘
                                                            ▼ Problem(case, solver)
                  ┌─ Problem: what and how, built ───────────────────────────────────────────────────────┐
                  │ Problem(case, solver)                    always both: Problem(case, Solver(case))    │
                  │                                                                                      │
                  │ mpts   ::Point   from case_physics (zones, constants, laws, grf) and solver.dtype    │
                  │ mesh   ::Mesh    from the case (L, nel, bcs)                                         │
                  │ basis  ::Basis   from mpts, mesh and solver.basis                                    │
                  │ time   ::Time    from solver.run (te, tg, tep)                                       │
                  │ solver ::Solver  as given                                                            │
                  └──────────────────────────────────────────────────────────────────────────────────────┘
                                                            ▼
                  ┌─ Simulation: a problem ready to run ─────────────────────────────────────────────────┐
                  │ Simulation(problem; folder)                                                          │
                  │ setup(case, solver; folder): builds one, saves it                         │
                  │                                                                                      │
                  │ problem ::Problem                                                                    │
                  │ folder  ::String   where its files go; by default the case's                         │
                  └──────────────────────────────────────────────────────────────────────────────────────┘
                    execute!(sim) runs solver.run.workflow; a .jld2 is a saved Simulation
```

Until `refactor-case-in-easydraw` ElastoPlasm had its own `Case` (the layout and settings with their typed reading as
fields, `case.material`, `case.bcs`, …, and a keyword constructor `Case(; L, nel, …)`); the user asked for "Case,
Layout and Setting [to] belong to EasyDraw, whereas Problem is from ElastoPlasm", towards one goal: "make things
independent of each other wherever possible". Until `refactor-solver-split` the solver held its case (`solver.case`) and copies of the case's laws, random field and
boundaries, and `Problem(solver)` built from it; the user chose the other direction ("Problem unifies what and how").
Every type's docstring lists its fields (`# Fields`), in the house layout of `conventions.md`.

- `Mesh{T1,T2,D}` — Eulerian background grid. Carries no connectivity (`e2n`/`e2e` live
  on `Basis`).
- **Host or device arrays: the containers' array types are type parameters.** `Point`, `PointSolidPhase`, `Mesh`,
  `MeshSolidPhase`, the thermal phases, `Basis` and `ApicTransfer` take one parameter per distinct array type (`V`, `VV`,
  `VM`, …), appended after the parameters code dispatches on, so `Point{T1,T2,D,CM,ST,EL}`, `Mesh{T1,T2,D}`,
  `Basis{T1,T2,D,NN,K,TR}` still match every version and field syntax is unchanged. These parameters are unbounded:
  bounds such as `V<:AbstractVector{T2}` refer to the other parameters and made subtype checks, hence compilation,
  about three times slower (a first problem built 10 → 34 s, `test_case` 43 s → 9 min); partial constructors
  (`Mesh{T1,T2,D}(…)`, `Basis{T1,T2,D,NN}(…)`, …) take the arrays' types from the arguments. Each container's
  `Adapt.adapt_structure` rebuilds it through `adapt_fields` (`boot/needs/utils.jl`), keeping plain-bits fields, so
  `adapt(backend, mpts)` holds that backend's arrays and `adapt(Array, ·)` brings them back with the same type.
  `MeshProperties{T1,T2,D,N,M}` is a plain value (static vectors), read alike on the host and in kernels; `mpts.f`/`t`
  and `mesh.t` are `Nothing` or the phase through parameters, not `Union` fields; `basis.e2e` is an `Adjacency`
  (`ptr`/`idx`, `nzrange`/`rowvals` like a sparse matrix's). The CPU code is unchanged (bit-identical, single thread).
- `Point{T1,T2,D,CM<:AbstractConstitutiveModel,ST<:AbstractStrain,EL<:AbstractElasticLaw,L,…}`
  — material points (Lagrangian). Every parameter up to `EL` is a real dispatch axis; `L == D*D`
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
- `Solver{T1,T2,D,…}` (`src/boot/needs/types/solver.jl`) — how a case is solved: one struct, its `solution` the
  option's text, `"explicit"` or `"implicit"` (checked in `build_solver` against `get_option()`), which the explicit
  workflows check (`explicit_only`: an `ArgumentError` for an implicit solver). Fields: the numerical option sections, `run` (`solver_run`: te, tg, tep, workflow, lstps; no folder) and the kernels (`cairn`).
  `Solver(case; run=(;), overrides…)` (`src/home/init/case/simulation.jl`) reads a case's settings (the option sections
  that are not `get_physics_sections()`, the run settings), overridable per field, and keeps nothing of the case;
  `Solver(; dim, sections…)` (`src/home/api/solver/build.jl`) is one written in code, `get_default()` with the
  sections merged. Both refuse physics (`material`, `grf`, `bcs`: `refuse_physics`) and end with `build_solver`.
  `solution` follows the workflow unless set; an explicit workflow on an implicit solver is refused.
- `Problem{T1,T2,D,P,M,B,S}` (`src/boot/needs/types/problem/problem.jl`) — what is solved and how, built: `mpts`,
  `mesh`, `basis`, `time`, `solver`, built by `Problem(case, solver)`: the mesh from the case's domain and `bcs`
  (`setup_mesh(L, nel, bcs, solver)`), the points from its zones, constants, `laws` and `grf` (`case_points(case,
  mesh)`, `setup_mpts(mesh, pts, laws)`), in the solver's precision, the time stages from `solver.run`, the basis from
  `solver.basis`. A solver of another dimension is refused. It is the one constructor (the user's decision: "A Problem
  constructor always need a Case and a Solver"); a case solved as its settings say is `Problem(case, Solver(case))`.
  The parts' types are unbounded parameters, set by `Problem(case, solver)` from what it builds, all of the solver's
  `T1`, `T2`, `D` (no bounds naming other parameters, as in every container); `save_simulation`/`load_simulation` swap
  the solver of a built problem through `with_solver`. What workflows run on (`elastodynamic!(problem)`).
- `Simulation{P}` (`src/boot/needs/types/problem/simulation.jl`) — a problem ready to run in a folder: `problem`,
  `folder` (absolute; a name given is under `dump/`), built by `Simulation(problem; folder)`; its workflow and load
  steps are its solver's (`problem.solver.run`). `setup(case, solver; folder=case_folder(case))` builds,
  saves and returns one's `.jld2` (`setup(case)` with the case's own solver), and writes beside it the case
  as solved (`save_case(file, case, solver)`: the solver's numerical options that differ from `get_default()` and its
  run settings, so `Solver(load_case(file))` is the same solver). A `.jld2` is a saved one
  (`save_simulation`/`load_simulation`); `execute!(sim)` / `execute!(jld2)` run it in place (`execute` on a copy).
- `Case` (EasyDraw's, `lib/EasyDraw/src/case.jl`, imported and re-exported) — what is solved: a `Layout`, the
  geometry (`case.layout`: domain, grid, zones, each a name and a shape `Polygon`/`Block`/`Ball`/`Region`), and its
  `Settings`, the values (`case.settings`: zones' values by name, defaults, sides, options, run), its own copies,
  checked against each other (`check_settings`); EasyDraw gives the values no meaning. Its file is
  `save_case`/`load_case` (EasyDraw's; geometry and settings may also be separate files, `load_case(geometry,
  settings)`; a file with top-level tables a case file does not have, such as ElastoPlasm's older `bcs`/`material`/
  `solver`, is refused). **ElastoPlasm reads it** (`src/boot/needs/types/problem/case.jl`): `case_physics(case)`, the
  typed reading of its physics, checked (`material`: `get_material()` with the defaults; `properties`: each zone's
  constants by name, against `get_zone_properties()`; `bcs`: the sides, `side_axes`; `laws`: `(; elastic, plastic)`, the
  options' `material` section over `get_default()`'s; `grf`, likewise; `thermal`; it checks the numerical options and
  run settings too, which `Solver(case)` reads); `case_settings(; properties, material, bcs, solver, thermal, run)`,
  settings written in ElastoPlasm's terms (`Case(Layout(…), case_settings(…))`, what `Case(; …)` was); `case_folder(case)`
  (the run settings' `folder`, else the layout's file's folder, else `get_run().folder`); `case_points`, its
  `PointData`, the explicit per-point contract `setup_mpts` reads; `save_case(file, case, solver)`, a method of
  EasyDraw's. 2D or 3D (no extruded form: x–z shapes extend along y). The drawing window (EasyDraw's) edits it. See
  "Describing a case" and "Drawing a case" in `operations.md`.
- **The constitutive laws, by name and by type.** `get_constitutive_laws()` (`lagrangian.jl`, shaped like the `material` section) maps the option
  values (`"hencky"`, `"improved_hencky"`, `"hypoelastic"`; `"DP"`, `"VM"`) to the points' types (`ST`/`EL`; `CM`):
  `build_solid_phase`/`setup_cmp` read them one way (`table_value`), `laws(mpts)` the other, and `get_option`'s lists
  derive from them. After setup nothing reads the laws from the solver or the case: the run-time readers (the `σᵢ!`
  skip for hypoelastic in `mapsto`, plot file names, logs, HDF5 attributes, the initial plot, `sfrac.jl`) use
  `laws(mpts)`, and `init_mapsto` registers `σᵢ!` always.
- `Basis{T1,T2,D,NN,K<:AbstractBasisKind,TR<:AbstractTransfer}`
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

`ApicTransfer{T2,D,L}` carries `Bᵢⱼ`/`Dᵢⱼ` (`Vector{SMatrix{D,D,T2,L}}` — `L` is a separate
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
  `elastoquasistatic!` call the same `ignite(mpts,mesh,basis,solver::Solver)`
  wrapper, taking any solution (unlike `elastodynamic!`/`elastoplastic!`, which refuse an implicit
  solver, `explicit_only` — see `solution` in `operations.md`).
- **Running on a device (explicit path).** `elastodynamic!`/`elastoplastic!` move the points, mesh and basis to the
  solver's device once (`to_device`, `src/home/core/device.jl`: the same objects on the host, no copy; otherwise
  `adapt(exec_backend(solver), ·)`), step there, and copy the arrays back into the host problem (`to_host!`) at each
  checkpoint that plots or records and at the end: the host problem stays what a `Simulation` holds, saves and plots.
  A step touches the device arrays only through kernels and whole-array operations (`mapreduce`, `map!`, `fill!`,
  `copyto!`): `get_dt`'s wave speed and `mapsto`'s `vmax` are reductions, the nonlocal element buckets are built by
  kernels (`element_particles`). `elastoquasistatic!` and `thermodynamic!` refuse a device (`host_only`).
- `src/home/core/solver/explicit/` — the primary, actively-used path
  (`elastodynamic!`, `elastoplastic!`). Kernel dispatch table lives on
  `solver.cairn` (a `NamedTuple`, dot-accessed: `solver.cairn.ignite.tplgy!`,
  `solver.cairn.mapsto.map.p2n!`, `solver.cairn.update.deform!`, ...), built by
  `init_ignite`/`init_mapsto`/`init_update`/`init_implicit` in `api/solver/build.jl`, on the solver's backend
  (`exec_backend(instr)`, the device `backend.exec` selected; `CPU()` for the host); each launch is followed by
  `sync(bk)`, `bk = exec_backend(solver)`: the solver owns the backend, and the kernels, the syncs and `to_device`
  all take it from there (not from the arrays, `get_backend(mpts.x)`, as at first). Every
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
