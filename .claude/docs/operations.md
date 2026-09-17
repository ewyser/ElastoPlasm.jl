## Operating the package

Typical end-to-end run (this is also the standard smoke test after any refactor):

```julia
using ElastoPlasm
L,nel = [64.1584, 64.1584/4.0], [40,10]      # domain size, element counts (2D here; 3 entries = 3D)
jld2  = slump_problem(L, nel; cli()...)           # builds mesh/mpts/basis/time, saves setup, returns the .jld2 path
out   = elastoplasm!(jld2; workflows = [elastodynamic!, elastoplastic!])
@assert out.success
```

Use `elastoplasm!` (not `elastoplasm`) whenever you need to inspect post-run
`mpts`/`mesh` state afterward — see
`.claude/bug/fixed/mlsmpm-shpfun-dispatch-ambiguity.md` for why.

- `slump_problem(L, nel; fid="...", kwargs...)` (`src/home/script/example/slump.jl`) is
  the reference example problem. Calls, in order: `get_solver` → `setup_geometry` →
  `setup_problem` (internally: `setup_mesh` → `setup_material_constants` →
  `setup_mpts`, which calls `setup_cmp` to build `mpts.s.cmp` → `setup_time`; returns a
  `MechanicalProblem`) → `setup_basis`, then `export_problem(...)` to persist to
  `.jld2` and returns that path. Use it as the template for a new example problem.
- `collapse_problem(nel; kwargs...)` (`src/home/script/example/collapse.jl`) is a dry
  granular Drucker-Prager block-collapse example — a plain rectangular block released
  under gravity with no pre-existing slope, geometry via `get_collapse`
  (`src/home/init/mpts/get_collapse.jl`, masks the shared `mpts_populate` candidate
  grid down to a `[0,w0]×[0,h0]` box). Modeled directly on MaterialPointSolver.jl's
  `2d_collapse.jl` reference scenario (`LandslideSIM/Archive_MaterialPointSolver.jl_paper`)
  — domain/block size/material constants are hardcoded to match it, so `nel` is the
  only argument meant to normally change; everything else is an override. Uses the
  manual low-level `setup_geometry → setup_mesh → setup_material_constants →
  setup_mpts → setup_time → MechanicalProblem → setup_basis` pipeline (not the
  `setup_problem` wrapper), since it needs custom `ρ0`/`E`/`ν`/`ϕ`/`c0` that
  `setup_problem`'s fixed `(mesh,mat,solver)` geometry-helper signature can't pass
  through. Defaults `plast.status=true`/`plast.constitutive="DP"`, and — unlike the
  package default (`:roller`, frictionless normal-only slip) — a `:fixed` (sticky/
  no-slip) base boundary, needed for the material to actually pile up instead of
  sliding indefinitely; matches MaterialPointSolver.jl's own boundary treatment.
  **Not to be confused with `column_problem`** (`src/home/script/example/column.jl`,
  renamed from the old `collapse_problem`) — a 1-D elastic self-weight column
  convergence test with no plasticity at all; see
  `.claude/bug/fixed/test-collapse-renamed-and-fixed.md` for that history.
- `cli()` (`src/home/api/solver/cli.jl`) parses interactive/CLI overrides for solver
  config; `cli(; ui=true)` prompts interactively, `cli()` with no args picks defaults
  non-interactively — pass its result as `kwargs...` into `slump_problem`/`get_solver`.
- `get_solver(; dim=2, kwargs...)` (`src/home/api/solver/get_solver.jl`) merges
  user kwargs over `get_default()`, resolves the execution backend, and builds the
  `Cairn` dispatch table (`cairn.ignite`/`.mapsto`/`.update`/`.implicit`).
- `elastoplasm(jld2_path; workflows=[...])` reopens the `.jld2`, unpacks
  `(mesh, mpts, basis, time, solver)`, and runs each workflow function in order as
  `workflow!(mpts, mesh, basis, time, solver)`. `elastoplasm!` is the in-place variant
  (opens `"r+"` so postprocessing can write back). Both accept any callable matching
  that signature — the built-ins are `elastodynamic!`/`elastoplastic!` (explicit
  solver) and `elastoquasistatic!` (dynamic_relaxation), freely composable in the
  `workflows` vector since they share the same call signature.
- Plotting is opt-in via `solver.plot.status`; when on, `elastoplasm` auto-saves a PNG
  per workflow named `"$(dim)_$(basis.which)_$(solvertype)_$(deform)_$(workflow)_$(quantity).png"`
  under the run's `dump/.../plot` path.

### Running tests / benchmarks

```julia
julia --project=. -e 'using Pkg; Pkg.test()'          # or: include("test/runtests.jl")
```

`test/runtests.jl` interactively lets you pick which `test/testset/test_*.jl` files to
run (auto-runs all of them under `GITHUB_ACTIONS=true`). It swallows exceptions inside
`runtests()`'s `try/catch` (just increments a fail counter) — when debugging a failing
test file, `include()` it directly instead of through `runtests()` to see the real
stacktrace, e.g.:

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
partition-of-unity `Σᵢ Nᵢ(x)` — useful as a quick correctness check on a new/changed
basis kind. Both this file and `test_performance.jl` report numeric checks
informationally (`println`, not hard `@test` gates), deliberately — not every basis
kind is expected to hold exact PoU everywhere, so a strict `@test` would be a false
failure rather than a real regression signal. When adding LaTeX to a generated figure:
wrap genuine math expressions but keep plain identifiers (basis kind names) as plain
strings — a multi-letter plain word in math mode doesn't render reliably under the
default GR backend.

### Controlling solver behaviour via `defaults.jl`

`get_default()` (`src/home/api/solver/defaults.jl`) returns the base `NamedTuple`
config merged by `get_solver` — this is the single place that defines every tunable
knob and its out-of-the-box value.

- `solution` — `"explicit"`/`"implicit"`, picks the concrete solver struct at
  construction time (`ExplicitSolver`/`ImplicitSolver` — same field layout,
  `DynamicRelaxationSolver` unused scaffolding). Also what `elastoplasm.jl`'s
  generated filenames and `logs.jl`'s startup banner print. Picking `solution=
  "explicit"` then calling an implicit-only workflow fails with `MethodError` at the
  workflow call, not silently — `elastodynamic!`/`elastoplastic!` are typed to
  `ExplicitSolver` specifically. `ignite()` is the one place dispatching on the
  abstract `AbstractSolver` instead, since both paths call it regardless of `solution`.
- `dtype` — arithmetic precision (`bits`, element types `T0`)
- `basis` — `which` (`"bsmpm"`/`"gimpm"`/`"smpm"`/`"mlsmpm"`, see `get_basis`), `how`
  (GIMP domain update mode), `trsfr` (P2G/G2P transfer scheme: `"std"`/`"tpic"`/
  `"apic"`, see `get_transfer`), `C_pf` (PIC/FLIP blend). `trsfr`/`C_pf` live here
  rather than a separate `transfer` section since `Basis` (the struct) owns both
  `kind` and `transfer` as sibling fields — see "Transfer scheme dispatch" in
  `architecture.md`.
- `strain` — `deform` (`"finite"`/`"infinitesimal"`)
- `stab` — `locking` (F-bar volumetric locking correction on/off), `damping`, `musl`
  (MUSL velocity reprojection on/off — lives here rather than under `basis`/`transfer`
  since it's a stabilization technique applied regardless of transfer scheme).
- `bcs` — `dirichlet` boundary condition matrix, one `[lower upper]` row per dimension
- `grf` — Gaussian random field generator for heterogeneous cohesion/friction fields
  (`status` toggles it on; see `GRF.jl`)
- `plast` — `status`, `constitutive` (`"DP"`/`"VM"`/`"MC"`/`"camC"` — not all are
  wired up, check `setup_cmp`'s branch before relying on one; the `retmap` kernel
  dispatch lives on `Point`'s `CM` type parameter, not a runtime string — see "DP/J2
  retmap kernel unification" in `planned-improvements.md`)
- `nonloc` — non-local plastic strain regularization (`status`, `ls` length scale)
- `plot` — `status`, `freq` (plot every N `Time` checkpoints), `dpi`, `what` (list of
  field specs to plot, keyed by name via `get_mpts_variable_config()`)
- `perf` — `status`; when true, forces `deform="infinitesimal"` and disables `nonloc`
  and swaps in the `_fast` kernel variants for a lighter-weight run
- `backend` — `select` (execution backend, `"host"`/GPU target) and `distributed`

Two ways to change solver behaviour:

1. **Per-run override (preferred, no code change)** — pass any of these keys as
   keyword arguments through `slump_problem`/`get_solver`; `get_solver` only keeps
   kwargs whose keys already exist in `default` (unrecognized kwargs get a `@warn`,
   not an error), then merges them over `default` via `merge(ref, user)`. **This is a
   shallow-per-key merge** — to override one field of a nested `NamedTuple` block
   (e.g. just `basis.trsfr`) you must pass the *whole* section
   (`basis = (; which=..., how=..., trsfr=..., C_pf=...)`), not just the one field,
   since `merge` replaces the whole `:basis` entry rather than recursing into it —
   this has caused real, previously-live bugs (`thermal_problem`'s and
   `test_basis.jl`'s own hardcoded `basis=(;which=...)` overrides silently dropped
   `trsfr`/`C_pf` until fixed to spread `get_default().basis...` first). Example:
   ```julia
   slump_problem(L, nel; basis=(;which="gimpm",how="Uii",trsfr="apic",C_pf=1.0), strain=(;deform="finite"), stab=(;locking=true,damping=0.1,musl=true))
   ```
2. **Change the package-wide default** — edit the literal value in `get_default()`
   directly. Do this only for a genuine change of the shipped default behaviour, not
   for one-off experimentation (use kwargs for that).

New solver-config knobs must be added as a new key inside the relevant block (or a new
top-level block) in `get_default()`'s `default` NamedTuple, then threaded through
wherever `init_ignite`/`init_mapsto`/`init_update`/`init_implicit` branch on
`instr[:section][:key]` to select kernels — grep those `init_*` functions for the
existing pattern before adding a new one. If the new field belongs on the solver
struct itself (like `solution` did), it also needs threading through
`ExplicitSolver`'s/`ImplicitSolver`'s field list (`src/boot/needs/types/solver.jl`) and
the corresponding positional arg in `get_solver.jl`'s final
constructor call — both structs are kept in lockstep by hand.

### Using `cli()`

`cli(; ui::Bool=false)` turns `get_default()`'s config tree into a `Dict{Any,Any}` of
kwargs suitable for splatting into `slump_problem`/`get_solver` as `cli()...`.

- `cli()` (default, `ui=false`) — **non-interactive**, returns `get_default()`'s values
  reshaped into the kwargs `Dict`. What you want for scripts/automated runs.
- `cli(ui=true)` — **interactive**: first a `MultiSelectMenu` to pick which top-level
  sections to configure, then a `RadioMenu`/`MultiSelectMenu` per leaf option (nested
  sections handled recursively; sections with a `status` field skip sub-prompts when
  answered `false`; `plot.what` gets multi-selection). Unselected sections fall back
  to `get_default()`. Use only from a real terminal/REPL — blocks on
  `TerminalMenus.request`, never call from a non-interactive script or agent session.
- To override specific options while using `cli()` for everything else, splat `cli()`
  first then override afterward — later kwargs win. Same shallow-merge caveat applies.
- `get_option()` is a reference for valid values per key (e.g.
  `get_option().basis.trsfr` → `["std", "tpic", "apic"]`) — some tunables (e.g.
  `plast.constitutive`) accept values it doesn't fully enumerate; check `init_update`'s
  dispatch in `update/update.jl` if in doubt whether a value is actually wired up.
