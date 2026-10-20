## Operating the package

Typical end-to-end run (this is also the standard smoke test after any refactor):

```julia
using ElastoPlasm
L,nel = [64.1584, 64.1584/4.0], [40,10]      # domain size, element counts (2D here; 3 entries = 3D)
jld2  = slump_problem(L, nel; cli()...)           # builds the Problem, saves it, returns the .jld2 path
out   = elastoplasm!(jld2; workflows = [elastodynamic!, elastoplastic!])
@assert out.success
p     = load_problem(jld2)                        # the post-run Problem: p.mpts, p.mesh, ...
# or entirely in memory, no file:
p     = elastoplasm(Problem(slump_case(L, nel; cli()...)); workflows=[elastodynamic!, elastoplastic!])
```

Use `elastoplasm!` (not `elastoplasm`) whenever you need to inspect post-run
`mpts`/`mesh` state afterward — see
`.claude/bug/fixed/mlsmpm-shpfun-dispatch-ambiguity.md` for why.

- Every example is a **case constructor** plus a thin wrapper: `slump_case(L, nel; kwargs...)` returns
  a `Case` (see "Describing a case" below), `slump_problem(L, nel; kwargs...) =
  setup_simulation(slump_case(...))` builds and saves it and returns the `.jld2` path
  (`src/home/script/example/*.jl`). `kwargs` are solver options, passed over each example's own
  defaults (a section given replaces the example's section, then `get_solver` merges it per field over
  `get_default()`). Write a new example as a case constructor too.
- `slump_case`: two polygons, the slope (`z = 0.625 Lx − 1.25 x`, up to the highest sub-cell boundary
  below `lz = 12.8`, i.e. whole sub-cells, as the old `get_slump` kept them) and the base layer below
  `0.3 lz` (`phi = 7.5°`, three times warmer); 2D or 3D (extruded along y).
- `collapse_case(nel; kwargs...)`: a dry granular Drucker–Prager block `[0, w0] × [0, h0]` released
  under gravity, as MaterialPointSolver.jl's `2d_collapse.jl` reference
  (`LandslideSIM/Archive_MaterialPointSolver.jl_paper`): domain, block and material (`ρ0`, `E`, `ν`,
  `ϕ`, `c0 = cr`) match it, `nel` is what normally varies; `:fixed` (no-slip) base, needed for the
  material to pile up. **Not to be confused with `column_case`**, the elastic self-weight column
  convergence test (`test_column.jl`); see `.claude/bug/fixed/test-collapse-renamed-and-fixed.md`.
- `block_case(nel=[64,36]; kwargs...)` is the 8 × 8 m self-weight block of Sutcliffe et al. (2026)
  §4.3, a `collapse_case` with the paper's values as overridable defaults (von Mises `c0 = 20 kPa` Lode
  reading, smooth base, GIMP). Run it with `elastoquasistatic!(a...; lstps=10, G=10.0,
  plastic=false|true)`; elastic, the base pressure is 37.7 kPa (~16 s).
- `column_case(nel, ν, E, ρ0, l0)`: one `Block` up to `l0` in a `1.25 l0` domain one element wide,
  gravity ramped over `tg = ⌈40 · 2l0/c⌉`; 2D, or 3D with `nel = [1, 1, nz]` (impossible before `Case`).
  `thermal_case`: the whole domain, `thermal = true`, fixed sides. `collision_case` (`ic_collision`):
  two `Ball`s moving towards each other.
- The old `get_slump`/`get_collapse`/`get_column`/`get_thermal`/`get_collision` generators,
  `material_fields`, `setup_problem` and `setup_material_constants` are gone: a shape no polygon,
  block or ball describes is a `Region(x -> …)` zone, so there is one construction path.
- `cli()` (`src/home/api/solver/cli.jl`) parses interactive/CLI overrides for solver
  config; `cli(; ui=true)` prompts interactively, `cli()` with no args picks defaults
  non-interactively — pass its result as `kwargs...` into `slump_problem`/`get_solver`.
- `get_solver(; dim=2, kwargs...)` (`src/home/api/solver/get_solver.jl`) merges
  user kwargs over `get_default()`, resolves the execution backend, and builds the
  `Cairn` dispatch table (`cairn.ignite`/`.mapsto`/`.update`/`.implicit`).
- `elastoplasm(problem::Problem; workflows=[...])` runs each workflow in order, in memory, and
  returns the problem: `workflow!(problem)` when the workflow takes one (the built-ins
  `elastodynamic!`/`elastoplastic!`/`thermodynamic!` (explicit) and `elastoquasistatic!`
  (dynamic relaxation) do, and so do closures such as `(a...) -> elastoquasistatic!(a...;
  lstps=40)`), otherwise `workflow!(mpts, mesh, basis, time, solver)` (`run_workflow!`).
  `elastoplasm(jld2)` does the same on the problem saved in a file and discards the result;
  `elastoplasm!(jld2)` writes the problem back after each workflow (`load_problem` to inspect it).
  The workflows are freely composable in the `workflows` vector.
- Plotting is opt-in via `solver.plot.status`; when on, `elastoplasm` auto-saves a PNG
  per workflow named `"$(dim)_$(basis.which)_$(solvertype)_$(deform)_$(workflow)_$(quantity).png"`
  under the run's `dump/.../plot` path.

### Describing a case (`Case`)

A problem starts from a `Case` (`src/boot/needs/types/problem/case.jl`), plain data built by one path,
`setup_simulation` (`src/home/init/case/`), stored as TOML.

```julia
case = Case(; L=[20, 10], nel=[20, 10], ni=2, bcs=[:roller :roller; :fixed :free],
              zones=[Zone("soil", Block([0, 0], [20, 6])),
                     Zone("rock", Polygon([[0, 0], [20, 0], [20, 2], [0, 2]]); E=5e7, nu=0.2, rho0=2900.0, phi=30*π/180)],
              material=(; c0=10.0e3), solver=(; basis=(; which="gimpm")), run=(; fid="layers", te=10.0, tep=5.0))
check_case(case)                      # errors, warnings, points/elements, estimated time step
jld2 = setup_simulation(case)         # dump/layers/layers_simulation.jld2, with layers_case.toml beside it
out  = run_simulation(case)           # setup + elastoplasm! with the run's workflow (case_workflows)
save_case("layers.toml", case); case = load_case("layers.toml")
```

- **Fields:** `L`, `nel` (2D or 3D; lower corner at the origin), `ni` points per element and axis,
  optional `width`/`nely` (a 2D x–z case extruded along y, y sides roller), `zones`, `material` (overrides
  of `MATERIAL`: E, nu, rho0, c0, cr, phi, Hp, n0, T, cp, k; SI, angles in radians), `bcs` (D×2
  `:free`/`:roller`/`:fixed`; extra rows ignored, as the solver does in 2D; `bcs` passed among the solver
  options moves here), `solver` (overrides, validated against `get_default()` at construction),
  `thermal`, `run` (`RUN`: `fid`, `te`, `tg`, `tep`, `workflow` ∈ `WORKFLOWS`, `lstps`).
- **Zones:** `Zone(name, shape; props...)`, props any `MATERIAL` key or `v` (initial velocity; two
  (x, z) components in 3D). Shapes: `Polygon` (x–z, even–odd, extruded in 3D), `Block(lo, hi)` (bounds
  included; D coordinates, or two for an extruded block), `Ball(centre, r)` (disk, sphere, or cylinder
  along y), `Region(f)` (a predicate on the point, code only, not writable to TOML). A point belongs to the **last** zone containing its sub-cell
  centre (`case_points`); unset constants come from `case.material`; the random field `solver.grf`
  perturbs the cohesion of points whose zone does not set `c0`, evaluated over all points at once (the
  RNG sequence of the old generators).
- **`PointData`** is the explicit per-point contract (`xp`, `vp`, `zone`, `c0`, `cr`, `phi`, `Hp`, `E`,
  `nu`, `rho0`, `n0`, `T`, `cp`, `k`) that `setup_mpts`/`build_solid_phase`/`build_thermal_phase`/
  `setup_cmp` read; `PointData(xp, ni)` is a uniform default material (used by `sfrac.jl`'s
  quadrature points). `mpts_populate(props; ni)` gives the candidate positions both use.
- **TOML** (`save_case`/`load_case`): `[domain]`, `[material]` (only values differing from `MATERIAL`;
  `phi` in degrees, rounded to 12 digits), `bcs`, `thermal`, `[solver.<section>]` (booleans, numbers and
  number lists as such, other values as their `get_option()` label, `plot.what` as variable names;
  `option_value` maps back and accepts any number for a numeric option), `[run]`, `[[zone]]` with
  `shape = "polygon"`/`"block"`/`"ball"`.
- **Checks** (`check_case`): errors (no material points, `solution` contradicting the workflow,
  quasi-static with hypoelastic), warnings (zone outside the domain, zone without points, elongated
  cells under dynamic relaxation, `tg > te`, 32-bit with GIMP), info (points, elements,
  `estimate_dt`: `get_dt` at rest, tested equal to it, and the number of steps).
- Verified: every example's setup bit-identical to the old generators (`test_setup.jl`, a stored
  digest reference; the collision as a point set, its disks now interleave); the standard slump run
  bit-identical; the §4.3 block under dynamic relaxation bit-identical single-threaded (with 8 threads it
  varies by ~1e-10 m between runs of the same code: atomic summation order).

### Recording a run (HDF5 output)

```julia
case = slump_case(L, nel; cli()..., plot=(; status=false),
                  output=(; status=true, every=0.5, what=(; point=(; solid=["u", "p", "epII"]))))
jld2 = setup_simulation(case)
elastoplasm!(jld2; workflows=[elastodynamic!, elastoplastic!])     # writes <name>_output.h5 beside the .jld2
f    = replace(jld2, "_simulation.jld2" => "_output.h5")
output_series(f, "point/solid/p"; point=30)                         # (; t, values, loadstep)
plot_output(f, "point/solid/p"; frame=7); animate_output(f, "point/solid/p"); plot_series(f, "point/solid/p"; point=30)
```

- **Configuration:** the `output` solver section, `(; status, every, single, what)`; `what` mirrors
  the data model, `point = (; base, solid, fluid, thermal)` and `mesh = (; solid, thermal)`, each a
  list of names from `OUTPUT_FIELDS` (`src/home/api/io/output/fields.jl`, a nested `Dict` shaped like
  `what`, one entry per quantity, extensible: point
  `x J n Omega ell`, solid `u v p dev epII epV rho q umag vmag`, thermal `T`; mesh solid `v a m`,
  thermal `T`). Partial selections merge per field; a phase the problem lacks is skipped with a
  warning. Multi-selections are any non-numeric option whose default is a list (`multi_select`):
  `plot.what` and these lists, checkbox grids in the window (pages "output · point/mesh").
- **File:** one array per field, frames last (`/point/<phase>/<field>`, `/mesh/<phase>/<field>`),
  chunked per frame, deflate 3, Float64 or Float32 (`single`); `/frames/t` (a time, or a load step:
  `/frames/loadstep` 1), `/frames/step`, `/frames/workflow` (the workflow's name); the grid,
  the initial point properties and attributes (case TOML, version, laws). Plain HDF5 (h5py etc.).
- **When:** explicit workflows stop on the `plot.freq` grid and the `output.every` grid
  (`checkpoints`, `explicit/get.jl`; without output exactly the old stops, so results are bit-identical with output
  off); `elastoquasistatic!` records each load step; `elastoplasm` records the initial state and the
  state each workflow ends on. The recorder lives in `task_local_storage` for the run (`record!`
  does nothing without one), so workflow signatures and closures are unchanged. `elastoplasm!`
  appends (it continues the saved state), `elastoplasm` and `elastoplasm(problem; output=file)` start
  a new file, `setup_simulation` removes a stale one. Adding stops changes the time-step sequence
  (`get_dt` clips at each stop), so a run with a finer `every` than `plot.freq` differs slightly
  from one without output.
- **After the run:** `output_series`, `output_frame`, `load_output` read it; `plot_output`,
  `animate_output`, `plot_series` draw from the file alone (live-plot colour maps, units and limits
  where the field matches a plot variable). Checked: a frame drawn from the file matches the live
  plot of the same time; `test_output.jl` covers shapes, cadence, last frame = saved state, append,
  Float32, skipped phases, load steps and the plotting.

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

`test_setup.jl` compares every example's initial problem (slump 2D/3D/random field/fine, collapse,
block, column, thermal, collision) with `test/dataset/setup_reference.jld2`: per array its size and a
64-bit FNV-1a digest of its bits, so bit-exact without storing the arrays. Run it after any `init/`
change; regenerate the reference (delete the file) only for an intended setup change, and say why.
`test_case.jl` covers `Case`: shapes, validation, zone membership and overrides, TOML round
trip, `check_case`/`estimate_dt` (= `get_dt`, 2D and extruded), zone constants reaching the particles,
free option values, workflows and a short `run_simulation`.

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
- `basis` — `which` (`"bsmpm"`/`"gimpm"`/`"smpm"`/`"mlsmpm"`, see `get_basis`), `trsfr` (P2G/G2P transfer scheme: `"std"`/`"tpic"`/
  `"apic"`, see `get_transfer`), `C_pf` (PIC/FLIP blend). `trsfr`/`C_pf` live here
  rather than a separate `transfer` section since `Basis` (the struct) owns both
  `kind` and `transfer` as sibling fields — see "Transfer scheme dispatch" in
  `architecture.md`.
- `stab` — `locking` (F-bar volumetric locking correction on/off), `damping`, `musl`
  (MUSL velocity reprojection on/off — lives here rather than under `basis`/`transfer`
  since it's a stabilization technique applied regardless of transfer scheme).
- `bcs` — `dirichlet` boundary condition matrix, one `[lower upper]` row per dimension
- `grf` — Gaussian random field for heterogeneous initial cohesion (`status` toggles it on;
  `covariance` `"gaussian"` or `"exponential"`; `param`: standard deviation `σ`, per-axis correlation lengths
  `Iₓ`, harmonics `Nₕ`, max wave number `kₘ`). Evaluated at the material points' own coordinates
  by `property_field`/`random_field` (`init/mpts/GRF.jl`), following Räss et al. (2019), in
  `refs/ElastoPlasm/grf/`.
- `material` — `elastic` (`"hencky"`/`"improved_hencky"`/`"hypoelastic"` — picks both
  the strain kinematics and the elastic law, see `AbstractElasticLaw` in `lagrangian.jl`;
  `"hypoelastic"` is the small-strain Jaumann-rate path, the two Hencky laws are finite
  strain) and `plastic` (`"DP"`/`"VM"` — the only two `setup_cmp` builds; any other value
  throws at setup time; the `retmap` kernel
  dispatch lives on `Point`'s `CM` type parameter, not a runtime string — see "DP/J2
  retmap kernel unification" in `planned-improvements.md`)
- `nonloc` — non-local plastic strain regularization (`status`, `ls` length scale)
- `plot` — `status`, `freq` (plot every N `Time` checkpoints), `dpi`, `what` (list of
  field specs to plot, keyed by name via `get_mpts_variable_config()`)
- `backend` — `select` (execution backend, `"host"`/GPU target) and `distributed`

Two ways to change solver behaviour:

1. **Per-run override (preferred, no code change)** — pass any of these keys as
   keyword arguments through `slump_problem`/`get_solver`; `get_solver` only keeps
   kwargs whose top-level keys exist in `default` (unrecognized ones get a `@warn`, not
   an error), then merges them over `default` field by field (`merge_config`): a partial
   section changes only the fields it names (`basis=(; trsfr="apic")` keeps `which` and
   `C_pf`), nested sections recurse (`grf=(; param=(; σ=1.0e4))`), and anything that is not
   a section (values, vectors such as `plot.what`) replaces the default. A nested field
   the defaults lack is an error naming the valid ones (`material=(; elastik=…)`); a
   `dtype` with only `precision` gets its `T0`/`bits`. (Until the merge became per field,
   a partial section silently dropped the other fields, a real bug twice:
   `thermal_problem`'s and `test_basis.jl`'s `basis=(;which=...)` lost `trsfr`/`C_pf`;
   full sections, as most callers still pass, merge exactly as before.) Example:
   ```julia
   slump_problem(L, nel; basis=(; which="gimpm", trsfr="apic"), stab=(; damping=0.1))
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
  first then override afterward — later kwargs win (a later section replaces the
  earlier one in the kwargs, then merges per field over the defaults).
- `get_option()` is a reference for valid values per key (e.g.
  `get_option().basis.trsfr` → `["std", "tpic", "apic"]`) — but not every listed value is
  wired up: `backend.select` offers `"cuda"`/`"rocm"` while every kernel is hard-coded to
  `CPU()` (see the GPU item in `planned-improvements.md`).
