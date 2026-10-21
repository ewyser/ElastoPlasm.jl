## Operating the package

Typical end-to-end run (this is also the standard smoke test after any refactor):

```julia
using ElastoPlasm
L,nel = [64.1584, 64.1584/4.0], [40,10]          # domain size, element counts (2D here; 3 entries = 3D)
jld2  = setup_simulation(slump_case(L, nel; cli()...))   # Case → Problem → Simulation, saved; returns the .jld2 path
execute(jld2)                                       # runs its workflows, written back after each
p     = load_simulation(jld2).problem              # the post-run Problem: p.mpts, p.mesh, ...
# or entirely in memory, no .jld2:
sim   = execute(Simulation(Problem(slump_case(L, nel; cli()...))); workflows=[elastodynamic!, elastoplastic!])
```

**Three types, each built from the one above:** a `Case` (the description) → `Problem(case)` (what is
solved: mpts, mesh, basis, time, solver) → `Simulation(problem; workflow, lstps, folder)` (how it is
run, by workflow name, and where its files go). A `.jld2` is a saved `Simulation`.

- Every example is a **case constructor**: `slump_case(L, nel; kwargs...)` returns a `Case` (see
  "Describing a case" below), and `setup_simulation(slump_case(...))` builds and saves it
  (`src/home/script/example/*.jl`; the `*_problem` wrappers are gone). `kwargs` are solver options, passed over each example's own
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
  `thermal_case`: the whole domain, `thermal = true`, fixed sides. `collision_case`:
  two `Ball`s moving towards each other.
- The old `get_slump`/`get_collapse`/`get_column`/`get_thermal`/`get_collision` generators,
  `material_fields`, `setup_problem` and `setup_material_constants` are gone: a shape no polygon,
  block or ball describes is a `Region(x -> …)` zone, so there is one construction path.
- `cli()` (`src/home/api/solver/cli.jl`) parses interactive/CLI overrides for solver
  config; `cli(; ui=true)` prompts interactively, `cli()` with no args picks defaults
  non-interactively — pass its result as `kwargs...` into `slump_case`/`get_solver`.
- `get_solver(; dim=2, kwargs...)` (`src/home/api/solver/get_solver.jl`) merges
  user kwargs over `get_default()`, resolves the execution backend, and builds the
  `Cairn` dispatch table (`cairn.ignite`/`.mapsto`/`.update`/`.implicit`).
- `execute(sim::Simulation; workflows)` (`src/home/core/execute.jl`) runs each workflow in order on
  `sim.problem`, in memory, and returns `sim`; `execute(jld2; workflows)` does it on the file, the
  simulation written back after each workflow, and returns the path. By default the simulation's own
  workflow (`case_workflows(sim)`: its `workflow` name, the quasi-static ones in `lstps` load steps). A
  workflow is a function of a `Problem`: the built-ins `elastodynamic!`/`elastoplastic!`/
  `thermodynamic!` (explicit), `elastoquasistatic!` (dynamic relaxation), and closures such as
  `(a...) -> elastoquasistatic!(a...; lstps=40)`; freely composable in `workflows`. It replaced
  `elastoplasm(problem)`, `elastoplasm(jld2)` (which discarded its result), `elastoplasm!(jld2)` and
  `run_simulation`.
- Plotting is opt-in via `solver.plot.status`; when on, `execute` saves a PNG per workflow named
  `"$(dim)_$(basis.which)_$(solvertype)_$(deform)_$(workflow)_$(quantity).png"` into the
  simulation's folder.

### Describing a case (`Case`)

A problem starts from a `Case` (`src/boot/needs/types/problem/case.jl`), plain data built by one path,
`setup_simulation` (`src/home/init/case/`), stored as TOML. No GUI package needed for any of it.

```julia
case = Case(; L=[20, 10], nel=[20, 10], ni=2, bcs=[:roller :roller; :fixed :free],
              zones=[Zone("soil", Block([0, 0], [20, 6])),
                     Zone("rock", Polygon([[0, 0], [20, 0], [20, 2], [0, 2]]); E=5e7, nu=0.2, rho0=2900.0, phi=30*π/180)],
              material=(; c0=10.0e3), solver=(; basis=(; which="gimpm")), run=(; folder="layers", te=10.0, tep=5.0))
check_case(case)                      # errors, warnings, points/elements, estimated time step
jld2 = setup_simulation(case)         # dump/layers/layers_simulation.jld2, with layers_case.toml beside it
execute(jld2)                         # its workflow (case_workflows)
save_case("layers.toml", case); case = load_case("layers.toml")
```

- **Fields:** `L`, `nel` (2D or 3D; lower corner at the origin), `ni` points per element and axis,
  optional `width`/`nely` (a 2D x–z case extruded along y, y sides roller), `zones`, `material` (overrides
  of `MATERIAL`: E, nu, rho0, c0, cr, phi, Hp, n0, T, cp, k; SI, angles in radians), `bcs` (D×2
  `:free`/`:roller`/`:fixed`; extra rows ignored, as the solver does in 2D; `bcs` passed among the solver
  options moves here), `solver` (overrides, validated against `get_default()` at construction),
  `thermal`, `run` (`RUN`: `folder`, `te`, `tg`, `tep`, `workflow` ∈ `WORKFLOWS`, `lstps`; `folder` is where
  `setup_simulation` writes the simulation and the case TOML, named after its last part: a name under `dump/`
  or a path; older cases' `fid`/`dir` are read into it).
- **Zones:** `Zone(name, shape; props...)`, props any `MATERIAL` key or `v` (initial velocity; two
  (x, z) components in 3D). Shapes: `Polygon` (x–z, even–odd, extruded in 3D), `Block(lo, hi)` (bounds
  included; D coordinates, or two for an extruded block), `Ball(centre, r)` (disk, sphere, or cylinder
  along y), `Region(f)` (a predicate on the point, code only, not writable to TOML). `Block` is not named
  `Box` because GLMakie exports `Box` (ambiguous when both are loaded). A point belongs to the **last** zone containing its sub-cell
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

### Drawing a case

```julia
using ElastoPlasm, QML
case = draw_case()                       # opens on the domain step, blocks until Create…; returns a Case
case = draw_case(slump_case([64.16, 16.04], [40, 10]))   # edit a case (opens on the zones step)
case = draw_case("layers.toml")          # or a case file
execute(setup_simulation(case))
```

- The window is the weak-dependency extension `ext/DrawQMLExt/` (Qt Quick). `src/` keeps `draw_case`
  (`script/draw/draw_case.jl`) and the stub `drawing_window(case; step)`, whose fallback errors with
  the install command; the extension adds the `drawing_window(::Case)` method (`qml_window`).
  - `editing.jl` is everything that needs no Qt, included by `DrawQMLExt.jl` and, on its own, by
    `test_draw.jl` (so it is tested where QML is not installed): `editor.jl`, the state and every edit as
    a plain Julia `Editor` (`press!`/`move!`/`release!`/`key!`, `set_domain!`, `set_option!`, `done!`,
    `save!`/`load!`, …); `helpers.jl`, GUI-independent pieces (zone constants `PROPS`, the preview points
    `drawn_points`, `snap_point`, `insert_on_edge!`, `delete_vertex!`, the option rows); `view.jl`, the
    JSON the window reads (`background_json`, `overlay_json`, `panel_json`). `DrawQMLExt.jl` keeps only
    the QML glue. QML counts from 0, Julia from 1: every index crosses through `qml_index`/`julia_index`
    (`view.jl`), and the QML file does no index arithmetic (it gets the tab, `panel.tab`, not the step).
  - **Cost per mouse move.** The canvas has two layers (`background_json`/`overlay_json`): the
    background (grid, domain, sides, the zones' fills and the points over them) repaints after an edit,
    a zoom or a resize; the overlay (outlines, vertices, the polygon being drawn, the cursor, the
    handle) at every move. A move repaints only the overlay (the panel and `check_case` refresh after
    an edit, on release, on keys and buttons), and the cursor coordinates are drawn on the canvas
    (`readout`), not in the status line. The preview points, the checks and the points' JSON (pixels,
    0.1 px) are cached in `editor.cache` until what they depend on changes (`preview_points`, `checks`,
    `points_json`); while vertices are dragged the preview keeps the last points and follows on release.
    Measured per move (`move` + `scene`, 40 moves, a slope polygon): 628 points 8 → 1.1 ms, 15,700
    points 1,466 → 2.2 ms, 62,800 points 19,557 → 5.2 ms. Most of the old cost was quadratic: the
    points' pixel list was built with `reduce(vcat, …)` over one vector per point. The Qt side, timed in
    the window offscreen (`Date.now()` around each layer's paint): with one layer, every move re-parsed
    and redrew all the points, ~70 ms at 62,800 points (15 ms Julia and the string to Qt, 19 ms
    `JSON.parse` of 0.8 MB, 35 ms drawing); with the overlay, a move costs ~1 ms at any size (<0.5 kB).
  - `DrawCase.qml` is only a view: two `Canvas` layers draw `Julia.background()`/`Julia.overlay()`, the panel reads
    `Julia.panel()`, both JSON strings (nested vectors don't convert through QML.jl, and property-map
    updates arrive asynchronously), and each widget calls one `@qmlfunction` (`DrawQMLExt.jl`). The
    functions are registered by `qml_window` for each window (`register_functions`), not in
    `__init__`: QML.jl forgets them when a window's application ends, so a second `draw_case` in the
    same session got a `Julia` object without them ("Property 'panel' … is not a function").
  - **Install** QML in the default environment with Qt pinned to the project's version, since the
    project's Plots/GR already loads `Qt6Base_jll` 6.8.2 and a QML built on Qt 6.10 then fails with a
    missing `_qt_version_tag_6_10` symbol: `] activate; add QML Qt6Base_jll@6.8.2
    Qt6Declarative_jll@6.8.2 Qt6ShaderTools_jll@6.8.2 Qt6Svg_jll@6.8.2 Qt6Wayland_jll@6.8.2` (the set
    verified; resolves QML 0.12.1, jlqml_jll 0.9.0). Bump the pins together with the project's Qt. Without
    the window, `draw_case`'s error (`drawing_window_help`) says whether QML is missing or only not
    loaded, explains the `_qt_version_tag` failure, and gives the pins with the Qt version the session
    actually loads (`session_qt`, from `Base.loaded_modules`), so it stays right when Plots/GR move Qt.
  - Offscreen (tests, CI): `QT_QPA_PLATFORM=offscreen QT_QUICK_BACKEND=software
    QT_QUICK_CONTROLS_STYLE=Fusion`; `qml_window(case; shots=(dir, steps))` saves a PNG of each listed
    step and closes (`grabToImage` on a child item: the window's `contentItem` gives an empty image).
- It edits 2D (x–z) problems whose zones are polygons, two-coordinate blocks (opened as their
  rectangle, same points) or two-coordinate balls (circles: `DrawnRegion` with `centre`/`radius` and no
  vertices, kept as a `Ball`, so the points are the same; New zone ▸ Circle draws one, centre then a point on it;
  its centre and radius handle drag, and the panel types its centre and radius); a `Region` or a 3D
  problem is refused. The preview uses each zone's own `inside` rule (`shape_of`), so it equals
  `case_points` (checked on the collision case). Zone properties it does not
  show (`v`, `Hp`, `n0`, `cp`, `k`), the problem's material and its thermal flag are carried through.
- Four tabs: **Domain** (drag the upper-right corner handle, snapped to a tenth of the order of
  magnitude, `nice_step`, or type `Lx`, `Lz` and the element size `h`, square cells, `nel_for`; a
  width along y for 3D, `nely` from `h`; points per element; one drop-down per side), **Zones** (the
  editor below; c₀, cᵣ [kPa], φ [°], T, E [MPa], ν, ρ₀; φ hidden when the options pick
  `material.plastic = "VM"`; "+ New zone ▾" opens a menu of shapes, `new_zone!(e, shape)`: Draw (a polygon,
  vertex by vertex), Rectangle (two opposite corners; every corner lying on another zone's edge splits
  it), Triangle (three corners, closes by itself), Circle (centre, then a point on it); the vertex last
  clicked or dragged (`picked`, outlined on the canvas) has one row, x and z typed exactly and "×"
  deleting it, `set_picked!`/`remove_picked!`, which move or delete it in every zone sharing it, through
  undo; the vertices are not listed), **Options** (one row per `get_option()` field: checkbox, name,
  drop-down; numeric fields, `free_numeric`, take any typed value plus presets; picking or typing a
  value ticks the box, unticking restores the default; multi-selections such as `plot.what` are
  checkbox grids; the `dtype` fields are the one `precision` choice and
  `solution` follows the workflow; `option_rows`) and **Run** (workflow, `te`/`tg`/`tep` or the load
  steps; `check_case` live; then **Create…** and **Cancel**).
- A picture to draw on (Domain step, "Picture…": a scanned section or a map, any image Qt reads), drawn
  under the grid at an adjustable opacity (`Backdrop`: file, pixel size reported by a hidden QML `Image`,
  lower-left corner and width in metres, the height from the aspect). The order is picture → mesh →
  draw: the picture opens across the domain's width; its scale is its width in metres, typed, or a scale
  bar ("Scale bar": click its two ends, type its length, "Apply": `apply_scale!`, a uniform scale about
  the picture's corner); setting the scale fits the domain to the picture (`fit_domain!`: from the origin
  to its far corner, at the current element size), which "Fit domain to picture" repeats after a crop.
  Its corner (`x0`, `z0`, 0 by default) moves it relative to the domain, e.g. to crop part of it. Checked
  offscreen on a hand-rasterised 20 × 10 m section with a 5 m scale bar: the bar, measured on the
  picture shown 10 m wide, gives width 20 m and a 20 × 10 m domain. A drawing aid only: not saved with
  the case (TOML).
- Only the tabs and the status line are always visible. "Start from a case file…" (top of the Domain
  step) loads a case TOML (`load!`; replaces the whole state and clears undo). **Create…** (end of the
  Run step) runs `check_case` first (`ready!`: errors block it, warnings need a second click), then opens
  Qt's folder dialog (choose, or make with its New Folder button): `create!` sets the run's folder
  (`run.folder`), writes `<name>_case.toml` there and closes with the case;
  `setup_simulation` then writes the simulation (`.jld2`, plots, outputs) into that folder instead of
  `dump/<name>/` (a plain name, as for cases made in code). There is no Save…: the case file is written by
  Create…. The dialogs start in the current directory (`qml_window(…; startdir=pwd())`), not in `dump/`,
  which is emptied whenever the package recompiles.
- Editing: clicks snap to existing vertices, then edges, then the sub-cell lattice `h/ni`
  (`snap_point`; particles sit at sub-cell centres, so never on a vertex), within a snap radius of 1.5%
  of the visible area (`tol`), so it follows the zoom (`zoom!`; the same radius picks vertices, edges
  and the domain corner). The cursor is grey when no zone is drawn, red while one is (with a crosshair
  pointer and "+ New zone" highlighted), and its shape says what it snapped to: a square on a vertex, a
  diamond on an edge, a ring on the lattice (`cursor_json`). While a zone is drawn, from
  its third vertex on, the cursor snaps onto its first vertex within the snap radius (a larger green
  ring, readout "click to close the zone"), and that click closes it, exactly there: no right click or
  Enter needed. With fewer vertices a click there is refused rather than repeating the vertex. The line's
  other vertices are not snap targets (the outline would touch itself). A vertex drawn onto an edge
  is also inserted into that edge, and dragging moves every zone's vertex at the same place, so shared
  boundaries stay shared (closing a zone inserts its edge-snapped vertices themselves,
  `insert_on_edge!(…; exact=true)`: projecting them again moved them by round-off, so the two zones'
  copies differed); clicking an edge inserts a vertex (`insert_on_edge!`), right-clicking deletes
  one (`delete_vertex!`); the wheel zooms about the pointer; undo/redo (Cmd/Ctrl+Z,
  Shift+Cmd/Ctrl+Z, buttons) keeps snapshots of the zones. The points to be created are previewed
  (`drawn_points`).

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
free option values, workflows, a short `execute`, and the `Simulation` (build, save, load, in memory equal to the file run, older files). `test_draw.jl` covers the drawing window: its Qt-free
part (`editing.jl`, loaded into a test module) everywhere, driven as a user would (drawing, snapping and
shared edges, drag, undo, domain, zone constants, options, run, Create… (its warning confirmation, the folder) and loading a case,
the caches and the readout, the JSON the view reads); then, where QML is installed, the window itself
offscreen, a PNG of each step, twice in one session.

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
  `DynamicRelaxationSolver` unused scaffolding). Also what `execute.jl`'s
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
   setup_simulation(slump_case(L, nel; basis=(; which="gimpm", trsfr="apic"), stab=(; damping=0.1)))
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
