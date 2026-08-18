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
- `Basis{T1,D,NN,K<:AbstractBasis}` (`src/boot/needs/types/concrete/basis/basis.jl`) —
  independent struct owning ALL mesh/point connectivity: `e2n`, `e2e`, `p2n`, `p2e`,
  `e2p`, `p2p`, plus `kind::K` (dispatches shape-function evaluation) and `NN` (nodes
  per element, only known once you fix mesh + basis kind, hence its own type param
  rather than living on `Mesh`/`Point`).
- Construction order is always **mesh → material points → basis** (basis depends on
  both). Function argument order is always **`(mpts, mesh, basis, ...)`** — this
  ordering is deliberate: it signals that the mpts↔mesh relationship is governed by
  basis.
- Shape-function evaluation is `eval_basis(mpts, mesh, basis, ip, nn)` (renamed from the
  old bare `basis(...)` to avoid colliding with the `Basis` type/struct field name).
  Concrete kinds live in `basis/{bsmpm,gimpm,smpm}.jl`; each defines
  `stencil_range(::ConcreteBasisKind, ::Type{T1})` used to build `Neighbs`.

## Explicit vs dynamic_relaxation solvers

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
  `ic_slump`/`elastoplasm(jld2; workflows=[...])` pipeline as the explicit solver.
  **Known gap**: `implicit.jl`'s `pt_solve_uP!` calls
  `instr.cairn.implicit.fint_p2n!(...)`, a key `init_implicit` never registers — the u-P
  path throws if actually exercised. Not fixed (out of scope until someone designs what
  `fint_p2n!` should do); `elastoquasistatic!` itself doesn't reach this code path.

**Known bug**: `basis.which="gimpm"` combined with `fwrk.locking=true` (F-bar volumetric
locking correction) crashes with `InexactError: trunc(Int64, NaN)`. Root cause:
`src/home/core/solver/explicit/update/fun/volumetric.jl`'s `ΔJp` divides unconditionally
by `mesh.s.m[no]` for every node `no` in a particle's basis stencil where the shape
function `N != 0`. GIMP's wider particle-domain support (`stencil_range` `-1:2`, reaching
further than bsmpm/smpm's compact 2-node support) means a particle's kernel can be
nonzero at a boundary/ghost node no particle has ever mapped mass to, so
`mesh.s.m[no] == 0` there and `mesh.ΔJ[no]/mesh.s.m[no]` evaluates to `0.0/0.0 = NaN`,
which propagates into `ΔFᵢⱼ` and corrupts particle positions. The analogous velocity
kernels (`mapsto/fun/solve.jl`, `mapsto/fun/augm.jl`) already guard this same division
pattern with an `iszero(mesh.s.m[no])` check; `volumetric.jl`'s `ΔJp` does not. Not fixed
(needs the same guard added, and needs deciding whether `ΔJn`'s corresponding
accumulation should also skip zero-mass nodes). Workaround: run gimpm with
`locking=false`, or use bsmpm/smpm when locking correction is needed.

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
  `"$(dim)d_$(solvertype)_$(deform)_$(workflow)_$(quantity).png"` under the run's
  `dump/.../plot` path.

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

`test_performance.jl` benchmarks core kernels directly via `solver.cairn.*` and compares
against `test/dataset/performance_baseline.jld2` (keyed by CPU name); delete/update that
file deliberately if a change is a genuine, intended perf shift rather than a regression.

### Controlling solver behaviour via `defaults.jl`

`get_default()` (`src/home/api/solver/defaults.jl`) returns `(solver_type, default)`
where `default` is the base `NamedTuple` config merged by `get_solver` — this is the
single place that defines every tunable knob and its out-of-the-box value:

- `dtype` — arithmetic precision (`bits`, element types `T0`)
- `basis` — `which` (`"bsmpm"`/`"gimpm"`/`"smpm"`, see `get_basis`), `how` (GIMP domain
  update mode: `"detFij"`/`"Fii"`/`"Uii"` finite, `"detΔFij"`/`"ΔFii"`/`"ΔUii"`
  infinitesimal), `ghost` (ghost-node count)
- `fwrk` — `deform` (`"finite"`/`"infinitesimal"`), `trsfr` (transfer scheme:
  `"std"`/`"tpic"`/`"apic"`), `C_pf` (PIC/FLIP blend), `musl` (MUSL velocity
  reprojection on/off), `locking` (F-bar volumetric locking correction on/off),
  `damping`
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
   (e.g. just `fwrk.trsfr`) you must pass the *whole* `fwrk = (; deform=..., trsfr=...,
   ...)` NamedTuple, not just the one field, since `merge` replaces the whole `:fwrk`
   entry rather than recursing into it. Example:
   ```julia
   ic_slump(L, nel; basis=(;which="gimpm",how="Uii",ghost=0), fwrk=(;deform="finite",trsfr="apic",C_pf=1.0,musl=true,locking=true,damping=0.1))
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
(`instr[:fwrk][:trsfr] == "apic"`-style branches) before adding a new one.

### Using `cli()`

`cli(; ui::Bool=false)` (`src/home/api/solver/cli.jl`) turns `get_default()`'s config
tree into a `Dict{Any,Any}` of kwargs suitable for splatting into `ic_slump`/
`get_solver` as `cli()...`. Two modes:

- `cli()` (default, `ui=false`) — **non-interactive**, just returns `get_default()`'s
  values reshaped into the kwargs `Dict`. This is what you want for scripts and
  automated runs (`ic_slump(L, nel; cli()...)`) — it's the plain-defaults path, not a
  prompt.
- `cli(ui=true)` — **interactive**: first shows a `MultiSelectMenu` letting you pick
  which top-level sections to configure (`dtype`, `basis`, `fwrk`, `grf`, `plast`,
  `nonloc`, `plot`, `perf`, `backend`); for each selected section it walks
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
  valid per key (e.g. `get_option().fwrk.trsfr` → `["std", "tpic", "apic"]`) — treat it
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
