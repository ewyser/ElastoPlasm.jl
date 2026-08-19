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
- `Point{T1,T2,D,E<:AbstractElasticity,R<:AbstractRheology,TM,TV,TS}` — material points
  (Lagrangian). No longer carries `NN`, `B`, or connectivity (`p2n`/`p2e`/`e2p`/`p2p`
  moved out). `mpts.x :: Vector{SVector{D,T2}}` — NOT a matrix; index fields with
  `getindex.(x, i)` or iterate, never `x[i, :]`.
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
  `ic_slump`/`elastoplasm(jld2; workflows=[...])` pipeline as the explicit solver — but
  is much slower to smoke-test (a full `elastoquasistatic!` run took ~5–8 minutes in
  practice, vs. seconds for `elastodynamic!`/`elastoplastic!`), so budget for that when
  verifying a change touches this path.

## Known bugs

Bugs discovered but not yet fixed, kept here so they don't get rediscovered from scratch.
**When picking one up: create a new branch off the current branch, named so the fix is
identifiable (e.g. `fix-volumetric-locking-zero-mass-guard`), rather than fixing in place
on an unrelated branch.**

- **`basis.which` ∈ {`"gimpm"`, `"mlsmpm"`} combined with `stab.locking=true` (F-bar
  volumetric locking correction) crashes with `InexactError: trunc(Int64, NaN)`.**
  Probable cause: `src/home/core/solver/explicit/update/fun/volumetric.jl`'s `ΔJp`
  divides unconditionally by `mesh.s.m[no]` for every node `no` in a particle's basis
  stencil where the shape function `N != 0`. Both GIMP and MLS use the wider `-1:2`
  particle-domain support (reaching further than bsmpm/smpm's node-type-corrected
  compact support) — GIMP because its kernel has no boundary correction at all, MLS
  because its moment-matrix projection is *designed* to replace node-type correction —
  so in both cases a particle's kernel can be nonzero at a boundary/ghost node no
  particle has ever mapped mass to. `mesh.s.m[no] == 0` there, so
  `mesh.ΔJ[no]/mesh.s.m[no]` evaluates to `0.0/0.0 = NaN`, which propagates into `ΔFᵢⱼ`
  and corrupts particle positions over subsequent timesteps (the crash itself surfaces
  later, in topology lookup, once positions go NaN — not at the division site). The
  analogous velocity kernels (`mapsto/fun/solve.jl`, `mapsto/fun/augm.jl`) already guard
  this same division pattern with an `iszero(mesh.s.m[no])` check; `volumetric.jl`'s
  `ΔJp` does not. Needs the same guard added, plus a decision on whether `ΔJn`'s
  corresponding accumulation should also skip zero-mass nodes. Workaround: run
  gimpm/mlsmpm with `locking=false`, or use bsmpm/smpm when locking correction is
  needed.
- **`dynamic_relaxation`'s u-P path throws if exercised.** Probable cause:
  `implicit.jl`'s `pt_solve_uP!` calls `instr.cairn.implicit.fint_p2n!(...)`, a key
  `init_implicit` never registers. Out of scope until someone designs what `fint_p2n!`
  should do; `elastoquasistatic!` itself doesn't reach this code path, so it isn't
  blocked by this.
- **`test/testset/test_collapse.jl` is fully non-functional**, beyond the
  `fwrk`→`strain`/`transfer`/`stab` kwarg fix already applied to it. Confirmed root
  cause by actually running it: **`ic_collapse` (`src/home/script/example/collapse.jl`)
  itself is broken** — it calls an undefined `kwargser` function and predates the
  current `get_solver`/`setup_geometry`/`setup_basis` pipeline entirely, the same
  vintage bug `ic_collision` had before this session's fix (see `collision.jl`'s
  rewrite for the pattern to follow: `get_solver` → `setup_geometry` → `setup_mesh` →
  `setup_cmpr` → `setup_mpts` → `setup_basis` → `setup_time` → `export_setup`).
  `ic_collapse` also never calls `setup_basis` at all (no `basis` in its `export_setup`
  call, which now requires one) and passes `misc` positionally to `export_setup`, which
  no longer accepts it (`export_setup` builds `misc` internally now). Once `ic_collapse`
  itself is fixed, three more bugs downstream in `test_collapse.jl` are already known
  and still need fixing: `load_simulation_setup` reads `file["cfg"]["instr"]`, the
  stale key name (`cfg["solver"]` is correct, see Persistence below);
  `compute_collapse_error` calls `elastoplasm!(jld2; workflow=[solver])` — singular
  `workflow`, not the real `workflows` kwarg; and `z0 = copy(setup.mpts.x[end, :])`
  matrix-style-indexes `mpts.x`, which is `Vector{SVector}` (see Core types).
- **`test/testset/test_slump.jl`'s full sweep (48 cases: 2 geom × 1 basis × 24
  `strain`/`transfer`/`stab` combos) fails 20/48, all under `locking=true`** —
  `infinitesimal`+`std`+`locking=true` cases in particular take 1+ minute each (vs
  ~15-20s for `locking=false` cases) before failing. Puzzling: a standalone
  single-case reproduction with the *identical* config (`bsmpm`, `infinitesimal`,
  `std`, `locking=true`) outside the test loop succeeded cleanly. Not yet root-caused
  — could be a real numerical issue specific to running many geometries/configs in
  sequence within one process (leftover global/mutable state between cases?), or
  something specific to the test harness's `@suppress`/`try`-`catch` interacting badly
  with a slow-but-eventually-successful run. Needs investigation with the full,
  untruncated test output (the run that surfaced this only had its last ~15 lines
  captured) before assuming it's the same class of issue as the already-diagnosed
  gimpm/mlsmpm F-bar locking NaN bug above (this was `bsmpm` only, which isn't
  supposed to hit that one).
- **`ic_collision`'s initial-velocity plot crashes**: `what_plot_field` errors with
  `type NamedTuple has no field data`, because `collision.jl`'s `plot.what` entries are
  hand-built as `(;mpts=(name="u",cblim=(...)))` instead of going through
  `get_mpts_variable_config()[name]` like `ic_slump`'s equivalent plot section does —
  the hand-built `NamedTuple` is missing the `data`/`label`/`unit`/`scale`/`cb` fields
  `what_plot_field` expects. `ic_collision` itself was also just rewritten this session
  (it previously called an undefined `kwargser` function and predated the current
  `get_solver`/`setup_geometry`/`setup_basis` pipeline entirely — that part is fixed),
  but this plot-config bug was pre-existing in the file and wasn't fixed. Workaround:
  run with `plot.status=false` until fixed.

## Planned improvements (not bugs, just known follow-up work)

- **Kernel granularity**: `elast!` is the model to follow — decomposed into small
  functions each doing one specific task (e.g. computing log strain), rather than one
  monolithic kernel body. Other large kernels (e.g. `update.jl`'s `elastoplast`/
  `elasto` dispatch functions, which branch on `strain.deform`/`basis.how` inline)
  would benefit from being split the same way.
- **Typed constitutive-model/tensor abstractions**, prototyped on branch
  `fancy-implementation` (`git show origin/fancy-implementation:<path>` — not merged,
  don't rebase onto it wholesale, see caveats below):
  - `lagrangian.jl`'s `AbstractConstitutiveModel{T,D}` with concrete `PerfectlyElastic`/
    `DruckerPrager` subtypes bundling *all* material parameters (elastic + plastic:
    `Gc`, `Kc`, `Del`, `Hp`, `c₀`, `cᵣ`, `ϕ₀`, `state`) into one typed object per
    particle (`cmp::Vector{CM}` on `PointSolidPhase`) — would replace the current split
    `E<:AbstractElasticity`/`R<:AbstractRheology` component design plus the loose
    `cmpr::NamedTuple` bag and `plast.constitutive::String` branch-dispatch in
    `update.jl`'s `init_update`.
  - `tensor.jl`'s `AbstractStrain`/`AbstractStress` typed tensors, each storing a
    **volumetric/deviatoric split** rather than a flat Voigt vector, with small
    dedicated functions per operation (`_trial_elastic_strain`, `_trial_elastic_stress`,
    `get_J2`, `get_τII`, `_get_cauchy_stress`, ...) dispatched on tensor *type*
    (`LogarithmicStrain`/`InfinitesimalStrain`, `CauchyStress`/`KirchhoffStress`)
    instead of on strings — a concrete instance of the kernel-granularity item above.

  Two caveats before adopting either: (1) that branch's `Point` embeds `basis::B`
  directly and uses a type-level `D<:AbstractDimension` instead of the current
  `D::Integer` — both predate/diverge from the current `Mesh`/`Point`/`Basis`
  separation and integer-dimension convention, so re-implement the constitutive-
  model/tensor ideas against the *current* architecture rather than merging that
  branch's `Point`/`Basis` shape backward; (2) the branch itself looks mid-refactor —
  `PerfectlyElastic` stores `Gc`/`Kc` as `Vector{T}` (one shared model, per-particle
  vectors) while `DruckerPrager` stores them as scalar `T` (one model instance per
  particle), an inconsistent convention for the same abstract type — so this isn't a
  drop-in copy-paste even for the parts worth keeping.
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
- **Check correctness and viability of the thermal solution.** `PointThermalPhase`
  has real fields (`c`, `k`, `q`, `T`) and several kernels reference it
  (`update/fun/heatflux.jl`, the `MeshThermalPhase` branches of `std_p2n`/`augm_p2n`/
  `n2p`), but those branches read `mpts.ϕ∂ϕ[nn,p,:]` — a field that does not exist
  anywhere on `Point` (confirmed by grep; see the caching-design note above for the
  fields that *do* exist). Any thermal-phase workflow would throw immediately if
  actually exercised — needs a real audit of whether this path was ever functional
  after the `Basis`/shape-function refactors, or needs rebuilding against the current
  `basis.N`/`∂N` cache the way the solid phase already was.

## Precision (32-bit) and StaticArrays performance gotchas

`dtype=(;T0=(Int32,Float32),bits=Int32(32),...)` is a supported but rarely-exercised
path — the default is always `(Int64,Float64)`, so bugs that only manifest under
`T1=Int32`/`T2=Float32` silently pass every day-to-day test. Two concrete patterns
found and fixed by actually running an `Int32`/`Float32` `ic_slump`→`elastoplasm`
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

Simulation setup is saved/loaded as `ic["mesh"]`, `ic["mpts"]`, `ic["basis"]`,
`ic["cmpr"]`, `ic["time"]`; solver config as `cfg["solver"]` (NOT `cfg["instr"]` — that
key name is stale/wrong despite appearing in some scripts). `elastoplasm`/`elastoplasm!`
unpack `mesh,mpts,basis,cmpr = file["ic/mesh"], file["ic/mpts"], file["ic/basis"],
file["ic/cmpr"]`. Both take a `workflows::Vector{Function}` kwarg (plural — not
`workflow`).

## Operating the package

Typical end-to-end run (this is also the standard smoke test after any refactor):

```julia
using ElastoPlasm
L,nel = [64.1584, 64.1584/4.0], [40,10]      # domain size, element counts (2D here; 3 entries = 3D)
jld2  = ic_slump(L, nel; cli()...)           # builds mesh/mpts/basis/cmpr/time, saves setup, returns the .jld2 path
out   = elastoplasm(jld2; workflows = [elastodynamic!, elastoplastic!])
@assert out.success
```

- `ic_slump(L, nel; fid="...", kwargs...)` (`src/home/script/example/slump.jl`) is the
  reference example problem. It calls, in order: `get_solver` → `setup_geometry` →
  `setup_mesh` → `setup_cmpr` → `setup_mpts` → `setup_basis` → `setup_time`, then
  `export_setup(...)` to persist everything to a `.jld2` file and returns that path.
  Use it as the template for defining a new example problem.
- `cli()` (`src/home/api/solver/cli.jl`) parses interactive/CLI overrides for solver
  config; `cli(; ui=true)` prompts interactively, `cli()` with no args picks defaults
  non-interactively — pass its result as `kwargs...` into `ic_slump`/`get_solver`.
- `get_solver(; dim=2, kwargs...)` (`src/home/api/solver/get_solver.jl`) merges
  user kwargs over `get_default()`, resolves the execution backend, and builds the
  `Cairn` dispatch table (`cairn.ignite`/`.mapsto`/`.update`/`.implicit`) — this is
  where solver-type kernels actually get instantiated (`SomeKernel(CPU())`).
- `elastoplasm(jld2_path; workflows=[...])` reopens the `.jld2`, unpacks
  `(mesh, mpts, basis, cmpr, time, solver)`, and runs each workflow function in order
  as `workflow!(mpts, mesh, basis, cmpr, time, solver)`. `elastoplasm!` is the in-place
  variant (opens the file `"r+"` so postprocessing can write back into it). Both accept
  any callable matching that signature as a "workflow" — the built-ins are
  `elastodynamic!`/`elastoplastic!` (explicit solver, `explicit/worflow.jl`) and
  `elastoquasistatic!` (dynamic_relaxation solver, `dynamic_relaxation/worflow.jl`);
  they can be freely composed/chained in the `workflows` vector since they all share
  the same `(mpts, mesh, basis, cmpr, time, solver)` call signature.
- Plotting is opt-in via `solver.plot.status` (set through `ic_slump`'s `plot=(...)`
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
   keyword arguments through `ic_slump`/`get_solver`; `get_solver` only keeps kwargs
   whose keys already exist in `default` (unrecognized kwargs just get a `@warn`, not
   an error), then deep-merges them over `default` via `merge(ref, user)`. Note this is
   a shallow-per-key merge — to override one field of a nested `NamedTuple` block
   (e.g. just `transfer.trsfr`) you must pass the *whole* `transfer = (; trsfr=...,
   C_pf=..., musl=...)` NamedTuple, not just the one field, since `merge` replaces the
   whole `:transfer` entry rather than recursing into it. Example:
   ```julia
   ic_slump(L, nel; basis=(;which="gimpm",how="Uii",ghost=0), strain=(;deform="finite"), transfer=(;trsfr="apic",C_pf=1.0,musl=true), stab=(;locking=true,damping=0.1))
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
tree into a `Dict{Any,Any}` of kwargs suitable for splatting into `ic_slump`/
`get_solver` as `cli()...`. Two modes:

- `cli()` (default, `ui=false`) — **non-interactive**, just returns `get_default()`'s
  values reshaped into the kwargs `Dict`. This is what you want for scripts and
  automated runs (`ic_slump(L, nel; cli()...)`) — it's the plain-defaults path, not a
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
  keyword-splat semantics: `ic_slump(L, nel; cli()..., basis=(;which="gimpm",
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
- Verify refactors by actually running the `ic_slump` → `elastoplasm(...; workflows=[...])`
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
