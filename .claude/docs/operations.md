## Operating the package

Typical end-to-end run (this is also the standard smoke test after any refactor):

```julia
using ElastoPlasm
L,nel = [64.1584, 64.1584/4.0], [40,10]          # domain size, element counts (2D here; 3 entries = 3D)
jld2  = setup_simulation(slump_case(L, nel; cli()...))   # Case + Solver → Problem → Simulation, saved; the .jld2 path
execute!(jld2)                                       # runs its workflows, written back after each
                                                     # (execute(jld2) runs a copy, returns it, leaves the file)
p     = load_simulation(jld2).problem              # the post-run Problem: p.mpts, p.mesh, ...
# or entirely in memory, no .jld2:
sim   = execute!(Simulation(Problem(slump_case(L, nel; cli()...))); workflows=[elastodynamic!, elastoplastic!])
# one case, another solver: its numerics overridden, the physics as it is
case  = slump_case(L, nel)
jld2  = setup_simulation(case, Solver(case; basis=(; which="gimpm"), run=(; te=5.0)); folder="slump_gimp")
```

**Four types:** a `Case` (what is solved: the geometry with its physics) and a `Solver` (how: numerical options, run
settings, kernels; `Solver(case)` reads a case's settings, overridable per field, `Solver(; dim, …)` is written in
code) meet in `Problem(case, solver)` (mpts, mesh, basis, time, solver), which `Simulation(problem; folder)` runs (where
its files go). The solver keeps nothing of the case: one solver solves several cases of its dimension, each written to
its own folder (`case_folder`). A `.jld2` is a saved `Simulation`. See the diagram in `architecture.md`.

- Every example is a **case constructor**: `slump_case(L, nel; kwargs...)` returns a `Case` (see
  "Describing a case" below), and `setup_simulation(slump_case(...))` builds and saves it
  (`src/home/script/example/*.jl`; the `*_problem` wrappers are gone). `kwargs` are solver options, passed over each example's own
  defaults (a section given replaces the example's section, then `solver_instr` merges it per field over
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
  reading, smooth base, GIMP). Its run workflow is `"elastoquasistatic!"` (so an implicit solver;
  `collapse_case` takes a `workflow` keyword): `execute!(jld2)` runs it elastic with `g = 9.81`; with the
  paper's gravity, `elastoquasistatic!(a...; lstps=10, G=10.0, plastic=false|true)`; elastic, the base
  pressure is 37.7 kPa (~16 s).
- `column_case(nel, ν, E, ρ0, l0)`: one `Block` up to `l0` in a `1.25 l0` domain one element wide,
  gravity ramped over `tg = ⌈40 · 2l0/c⌉`; 2D, or 3D with `nel = [1, 1, nz]` (impossible before `Case`).
  `thermal_case`: the whole domain, `thermal = true`, fixed sides. `collision_case`:
  two `Ball`s moving towards each other.
- The old `get_slump`/`get_collapse`/`get_column`/`get_thermal`/`get_collision` generators,
  `material_fields`, `setup_problem` and `setup_material_constants` are gone: a shape no polygon,
  block or ball describes is a `Region(x -> …)` zone, so there is one construction path.
- `cli()` (`src/home/api/solver/cli.jl`) parses interactive/CLI overrides for solver
  config; `cli(; ui=true)` prompts interactively, `cli()` with no args picks defaults
  non-interactively — pass its result as `kwargs...` into `slump_case` (not into `Solver(; …)`, which refuses its
  physics sections).
- `Solver(; dim=2, kwargs...)` (`src/home/api/solver/build.jl`) merges
  user kwargs over `get_default()` (`solver_instr`), resolves the execution backend, and builds the
  `Cairn` dispatch table (`cairn.ignite`/`.mapsto`/`.update`/`.implicit`), all in `build_solver`, which
  `Solver(case)` ends with too.
- `execute!(sim::Simulation; workflows)` (`src/home/core/execute.jl`) runs each workflow in order on
  `sim.problem`, in memory, and returns `sim`; `execute!(jld2; workflows)` does it on the file, the
  simulation written back after each workflow, and returns the path. `execute(sim)`/`execute(jld2)`
  run a copy and return it, leaving `sim` or the file as it was (`elastoplasm`'s logic; unlike
  `elastoplasm(jld2)`, the result is returned rather than discarded). By default the simulation's own
  workflow (`case_workflows(sim)`: its `workflow` name, the quasi-static ones in `lstps` load steps). A
  workflow is a function of a `Problem`: the built-ins `elastodynamic!`/`elastoplastic!`/
  `thermodynamic!` (explicit), `elastoquasistatic!` (dynamic relaxation), and closures such as
  `(a...) -> elastoquasistatic!(a...; lstps=40)`; freely composable in `workflows`. It replaced
  `elastoplasm(problem)`, `elastoplasm(jld2)`, `elastoplasm!(jld2)` and `run_simulation`.
- Plotting is opt-in via `solver.plot.status`; when on, `execute!` saves a PNG per workflow named
  `"$(dim)_$(basis.which)_$(solvertype)_$(deform)_$(workflow)_$(quantity).png"` into the
  simulation's folder.

### Describing a case (`Case`)

A problem starts from a `Case` (`src/boot/needs/types/problem/case.jl`), plain data, stored as TOML, solved by a
`Solver` (`setup_simulation`, `src/home/init/case/`). No GUI package needed for any of it.

```julia
case = Case(; L=[20, 10], nel=[20, 10], ni=2, bcs=[:roller :roller; :fixed :free],
              zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Polygon([[0, 0], [20, 0], [20, 2], [0, 2]]))],
              properties=Dict("rock" => (; E=5e7, nu=0.2, rho0=2900.0, phi=30*π/180)),
              material=(; c0=10.0e3), solver=(; basis=(; which="gimpm")), run=(; folder="layers", te=10.0, tep=5.0))
solver = Solver(case)                 # how: the settings' numerical options and run settings
check_case(case, solver)              # errors, warnings, points/elements, estimated time step
jld2 = setup_simulation(case, solver) # dump/layers/layers_simulation.jld2, with layers_case.toml beside it
execute!(jld2)                         # its workflow (case_workflows)
save_case("layers.toml", case); case = load_case("layers.toml")
```

- **A case is built from an EasyDraw `Layout` and its `Settings`** (branch `refactor-case-from-layout`, asked for: the
  case and the layout described the same model twice, then: "separate Layout (pure geometry) from Case, related to
  physics and solver"). The layout is the geometry: `L`, `nel` (2D x–z, or 3D x, y, z; lower corner at the origin),
  `ni` points per element and axis, `zones` (a name and a shape each, names unique). The settings are the values, as
  plain data: `zones` (each zone's values by zone name), `defaults`, `sides`, `options`, `run`. `Case(layout,
  settings)` checks a copy of both (`EasyDraw.check_settings`: values only for the layout's zones, sides its domain
  has) and holds the typed reading of its physics: `case.material`, `case.properties` (name → NamedTuple),
  `case.bcs`, `case.laws` (`(; elastic, plastic)`) and `case.grf` (both complete, over `get_default()`'s),
  `case.thermal`; the geometry is read from `case.layout`. The numerics (the other option sections, the run settings)
  are read by `Solver(case)` (branch `refactor-solver-split`: "Case is the combination between a geometry with physical
  properties assigned to it"), which keeps nothing of the case: `Problem(case, solver)` unifies them ("Problem
  unifies What and How"). A case keeps its own copies: to
  change a model, change a layout or settings and build a new case. One geometry takes several settings:
  `Case(layout, settings_a)`, `load_case(geometry_file, settings_file)`. `Layout`, `Settings`, `load_layout`,
  `load_settings`, `save_settings` are re-exported by `using ElastoPlasm`. `Case(; L, nel, ni, zones, properties,
  material, bcs, solver, thermal, run)` builds both, so the examples kept their shape (`slump` and `collision` moved
  their zone keywords into `properties`):
  - `properties`: by zone name, `Dict("base" => (; phi=…, T=…))`, any `MATERIAL` key or `v`; the settings' `zones`;
  - `material`: overrides of `MATERIAL` (E, nu, rho0, c0, cr, phi, Hp, n0, T, cp, k; SI, angles in radians), the
    settings' `defaults`;
  - `bcs`: D×2 `:free`/`:roller`/`:fixed` (extra rows ignored; `bcs` passed among the solver options moves here), the
    settings' `sides` through `side_axes(D)`: x `left`/`right`, z `bottom`/`top`, and in 3D y `front`/`back`, roller
    when a side is not set;
  - `solver`: option overrides, physics and numerics alike, the settings' `options` in their file form
    (`option_labels`; `Case` validates them all against `get_default()`, keeps the physics as `case.laws` and
    `case.grf`; `Solver(case)` reads the numerics). Every example's options round-trip exactly;
  - `thermal`, `run`: the settings' `run`, every `RUN` setting (`folder`, `te`, `tg`, `tep`, `workflow` ∈ `WORKFLOWS`,
    `lstps`) and `thermal` when true; `Solver(case)` reads all but the folder as `solver.run` (`solver_run`). The
    folder is the case's, `case_folder(case)`: where `setup_simulation` writes the simulation and the case TOML, named
    after its last part, a name under `dump/` or a path; without one in `run`, the layout's own folder (a file's). It
    is not the solver's, or one solver used on two cases would write both into one folder.
  - `Solver(case; run=(;), overrides…)`: numerical overrides field by field (a physics section, `bcs` or `run.folder`
    is refused), run settings with `run`. `setup_simulation(case, solver; folder)` writes beside the simulation the
    case as solved, `save_case(file, case, solver)`: its settings with the solver's numerical options (those differing
    from `get_default()`) and run settings, so `Solver(load_case(file))` is the same solver.
  - No `width`/`nely`: a 2D model extruded along y is a 3D case of x–z shapes (they extend along y), its y sides roller
    by default as the extruded form forced. EasyDraw has the same single 3D form.
- **Zones:** `Zone(name, shape)`, geometry only; a zone's constants are its `properties` (any `MATERIAL` key or `v`,
  initial velocity; two (x, z) components in 3D). `Zone`, `Layout`, `Settings` and the shapes are EasyDraw's
  (`lib/EasyDraw/src/{shapes,layout}.jl`, imported and re-exported in `boot.jl`), where a zone's values can be anything;
  `Case` checks them against `ZONE_PROPS` (an unknown one errors when the case is built, naming the zone). Shapes: `Polygon` (x–z, even–odd, extruded in 3D; several outlines, `Polygon(rings)`, make one zone of
  several parts or with holes: flat `vertices` plus `sizes`, `polygon_rings`), `Block(lo, hi)` (bounds
  included; D coordinates, or two for an extruded block), `Ball(centre, r)` (disk, sphere, or cylinder
  along y), `Region(f)` (a predicate on the point, code only, not writable to TOML). `Block` is not named
  `Box` because GLMakie exports `Box` (ambiguous when both are loaded). A point belongs to the **last** zone containing its sub-cell
  centre (`case_points`); unset constants come from `case.material`; the random field `case.grf`
  perturbs the cohesion of points whose zone does not set `c0`, evaluated over all points at once (the
  RNG sequence of the old generators).
- **`PointData`** is the explicit per-point contract (`xp`, `vp`, `zone`, `c0`, `cr`, `phi`, `Hp`, `E`,
  `nu`, `rho0`, `n0`, `T`, `cp`, `k`) that `setup_mpts`/`build_solid_phase`/`build_thermal_phase`/
  `setup_cmp` read; `PointData(xp, ni)` is a uniform default material (used by `sfrac.jl`'s
  quadrature points). `mpts_populate(props; ni)` gives the candidate positions both use.
- **TOML** (`save_case`/`load_case`): a case file is its layout's file with its settings (EasyDraw's
  `save_layout(file, layout, settings)`, `load_layout` + `load_settings`), SI throughout (φ in radians): the geometry
  at the top level, `[domain]` (`L`, `nel`, `ni`) and `[[zone]]` (name, `shape = "polygon"`/`"block"`/`"ball"`); the
  settings under `[settings]`: `[settings.zones.<name>]` (each zone's constants), `[settings.defaults]` (the material
  overrides), `[settings.sides]`, `[settings.options.<section>]` (booleans, numbers and number lists as such, other
  values as their `get_option()` label, quantities by key; `option_value` maps back and accepts any number for a numeric
  option), `[settings.run]` (with `thermal`). Geometry and settings can also be files of their own
  (`save_layout(file, case.layout)`, `save_settings(file, case.settings)`, `load_case(geometry, settings)`). A file of
  this branch's earlier format (values inline in `[[zone]]`, top-level `[sides]` …) still loads, through EasyDraw. `load_case` refuses a case file of the older format (`bcs`,
  `[material]`, `[solver]`, φ in degrees; the user chose dropping them over a converter): set it up again from its
  example. Every example rebuilds identically from its file (each array of the problem).
- **Checks:** `check_case(case)`, the geometry and physics: an error (no material points), warnings (zone outside
  the domain, zone without points), info (points, elements). `check_case(case, solver)` adds how it is solved: errors
  (`solution` contradicting the workflow, quasi-static with hypoelastic), warnings (elongated cells under dynamic
  relaxation, `tg > te`, 32-bit with GIMP), info (the load steps, or `estimate_dt`: `get_dt` at rest, tested equal to
  it, and the number of steps). An explicit workflow asked of an implicit solver is refused by `Solver(case)` itself.
- Verified: every example's setup bit-identical to the old generators (`test_setup.jl`, a stored
  digest reference; the collision as a point set, its disks now interleave); the standard slump run
  bit-identical; the §4.3 block under dynamic relaxation bit-identical single-threaded (with 8 threads it
  varies by ~1e-10 m between runs of the same code: atomic summation order).

### Drawing a case

```julia
using ElastoPlasm
case = draw_case()                       # opens on the domain step, blocks until Create…; returns a Case
case = draw_case(slump_case([64.16, 16.04], [40, 10]))   # edit a case (opens on the zones step)
case = draw_case("layers.toml")          # or a case file
execute!(setup_simulation(case))
```

- **The window is EasyDraw's** (`lib/EasyDraw/`, EasyDraw.jl, a package and repository of its own included as a git
  submodule, usable without ElastoPlasm, e.g. for finite differences or finite elements; its `README.md` and its documentation site, `docs/`, with a guide to
  the window and screenshots retaken by `docs/shots.jl`): it edits an `EasyDraw.Layout` (domain, zones, a boundary kind per
  side, the zone properties' defaults, options, run settings) with what a host's `EasyDraw.Catalog` offers (zone
  properties, boundary kinds, options, run settings as `Field`s: menu, numeric, integer, list; a display `scale`,
  `min`, `show_if`), TOML both. ElastoPlasm is one host: `get_catalog()` (`src/home/api/solver/catalog.jl`, beside `get_default`/`get_option`) builds its catalog
  from `get_option()`, `get_default()`, `MATERIAL` and `WORKFLOWS` (no hand-written copy) as plain data, a `Dict` shaped
  like the catalog's TOML, which EasyDraw reads (`EasyDraw.Catalog(dict)`; `draw` and `write_catalog` take either), so
  ElastoPlasm builds no EasyDraw type for it (the user's request: it first built `EasyDraw.Field`s, tying ElastoPlasm to
  EasyDraw's API rather than to the file format); EasyDraw stays a dependency for the shapes and `draw_case`.
  The window edits the case's layout and settings (`src/home/script/draw/draw_case.jl`): `draw_case(case)` is
  `Case(EasyDraw.draw(case.layout, case.settings; catalog=get_catalog(), case_hooks()...)...)`, no conversion.
  `case_hooks()` gives what a file cannot: the checks (`check_case(c, Solver(c))`, `c = Case(l, s)`), what Create… writes
  (`<name>_case.toml`, the chosen folder set as the run folder first: EasyDraw returns those same settings, so the case
  drawn has it too), what "Start from a file…" reads (a case file's layout and settings), the title.
  - **The window is not optional** (the user's decision: EasyDraw exists to draw interactively): QML is one of
    EasyDraw's dependencies, so `draw`/`draw_case` work without `using QML`. It was an extension (`EasyDrawQMLExt`, QML a
    weakdep the user installed in the default environment with Qt pins) until then. **QML loads with the first window, not
    with `using`** (`qml()` in `window/window.jl`: `Base.require(EasyDraw, :QML)`, its functions called through
    `invokelatest`): loaded with EasyDraw, it made `using ElastoPlasm` take 15–20 s instead of 4.6–5.8 s (QML itself
    ~8 s, mostly CxxWrap compiling, then ElastoPlasm's code recompiling, ~8 s more, as QML invalidates it); lazily,
    `using ElastoPlasm` is unchanged and `using EasyDraw` alone takes 0.01 s. A session's first window takes ~20 s
    (QML ~8 s, then compiling the window's code), the next ones 0.6 s; the very first in a new environment took ~2 min
    once (first-time caches). Checked in a fresh environment with EasyDraw alone (no ElastoPlasm): it resolves the same
    QML and Qt, and the window opens without `using QML`.
    In ElastoPlasm's environment the resolver keeps QML on Plots/GR's Qt (QML 0.12.1, jlqml_jll 0.9.0, every Qt6 library
    6.8.2), so the `_qt_version_tag` failure of a QML built for another Qt, which forced the pins while QML lived in
    another environment, cannot happen.
  - **Between the two, without loss** (asked for: "seamless compatibility"). *Points:* the window's preview,
    `EasyDraw.points`, `check_case`'s counts (`case_zones`) and the simulation's material points (`case_points`) are the
    same to the bit: one lattice, `EasyDraw.subcell_centres` (the range arithmetic of `mpts_populate`; the preview's own
    `(i - 0.5)h` formula differed by up to 7e-15 m, so a point on a zone's edge could change zone), and one membership
    rule, `EasyDraw.zone_at(zones, x)`, which `case_points` calls (setups bit-identical, `test_setup`). *Files:* a case
    file is a layout file, so `draw_case`, "Start from a file…" and `setup_simulation(load_case(…))` take a layout drawn
    with `get_catalog()` in EasyDraw alone, and EasyDraw opens a case file as it is. *Functions:* a `Case` works wherever
    EasyDraw takes a layout, through `case.layout` and `case.settings` (`zone_at`, also 3D, `zone_index`, `property` with the case's
    material as defaults, `points`, `partition`, `write_geo`, `save_layout`), exported by `using ElastoPlasm`; on a 3D
    case the x–z ones work in its x–z plane, while `check_case` counts the 3D points. *Option values:* a case file (`option_value`) and the
    window (`EasyDraw.field_value`) parse them separately, so `test_draw` checks they agree on every option over a set of
    inputs (438 cases); it found `option_value` accepting `Inf`/`NaN`, now refused as the window does.
  - Until `feat-easydraw-package` the window was ElastoPlasm's extension `ext/DrawQMLExt/`, editing a `Case` with
    ElastoPlasm's tables built in.
  - EasyDraw's `src/editor/` is everything that needs no Qt (tested without a window): `editor.jl`, the
    state and every edit as a plain Julia `Editor` (`press!`/`move!`/`release!`/`key!`, `set_domain!`, `set_option!`,
    `set_run!`/`pick_run!`, `ready!`/`create!`, `load!`, …); `helpers.jl`, GUI-independent pieces (the preview points
    `drawn_points`, `snap_point`, `insert_on_edge!`, `delete_vertex!`, the options' nested storage); `view.jl`, the
    JSON the window reads (`background_json`, `overlay_json`, `panel_json`). `src/window/window.jl` (beside
    `EasyDraw.qml`) keeps only the QML glue: one table, `CALLS`, of what QML may call by
    name, mostly `act(f)` (the open editor, QML strings as Julia strings, returns nothing), so a call cannot
    be defined and left unregistered; any parsing or logic lives in `editor.jl`, tested without Qt. QML counts from 0, Julia from 1: every index crosses through `qml_index`/`julia_index`
    (`view.jl`), and the QML file does no index arithmetic: it gets the tabs and the current one (`panel.tabs`,
    `panel.tab`; `tab_steps`: no Options tab when the catalog has no options) and the page (`panel.page`), and
    `set_step` takes a tab, which Julia maps to its step.
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
  - `EasyDraw.qml` is only a view: two `Canvas` layers draw `Julia.background()`/`Julia.overlay()`, the panel reads
    `Julia.panel()`, both JSON strings (nested vectors don't convert through QML.jl, and property-map
    updates arrive asynchronously), and each widget calls one of `CALLS` through `win.act(() => Julia.…)`,
    which refreshes the panel and canvas after it. The calls are registered (`QML.qmlfunction`) by
    `EasyDraw.window` for each window (`register_functions`), not in
    `__init__`: QML.jl forgets them when a window's application ends, so a second `draw_case` in the
    same session got a `Julia` object without them ("Property 'panel' … is not a function").
  - Offscreen (tests, CI): `QT_QPA_PLATFORM=offscreen QT_QUICK_BACKEND=software
    QT_QUICK_CONTROLS_STYLE=Fusion`; `draw(layout; shots=(dir, steps))` (or `draw_case(case; shots)`,
    `EasyDraw.window(editor; shots)`) saves a PNG of each listed step and closes (`grabToImage` on a child item: the window's `contentItem` gives an empty image).
    **Not what a user sees on macOS**: without the style override the window uses Qt's native macOS style,
    which follows dark mode and cannot render offscreen (it crashes in `QMacStyle`). Check a visual change
    there too, with `QT_QPA_PLATFORM=cocoa` and no style or backend override (the window shows for the shots
    and closes). Found that way: a fixed grey (`#444`) status line invisible on the dark panel, and a `Frame`
    without margins in that style (the picture's panel is a `GroupBox`); help texts are the palette's text
    at `opacity: 0.6`, not a fixed `gray`.
- It edits 2D (x–z) problems whose zones are polygons, two-coordinate blocks (opened as their
  rectangle, same points) or two-coordinate balls (circles: `DrawnRegion` with `centre`/`radius` and no
  vertices, kept as a `Ball`, so the points are the same; New zone ▸ Circle draws one, centre then a point on it;
  its centre and radius handle drag, and the panel types its centre and radius); a `Region` or a 3D
  problem is refused. The preview uses each zone's own `inside` rule (`shape_of`), so it equals
  `case_points` (checked on the collision case). Zone properties it does not
  show (`v`, `Hp`, `n0`, `cp`, `k`), the problem's material and its thermal flag are carried through.
- Four tabs (three without options): **Domain** (first a choice, *Domain* or *Picture to draw on*, the latter with its own panel,
  below; drag the upper-right corner handle, snapped to a tenth of the order of
  magnitude, `nice_step`, or type `Lx`, `Lz` and the cell size `h`, square cells, `nel_for`; a
  width along y for 3D, `nely` from `h`; points per cell; one drop-down per side, its kinds and colours the
  catalog's), **Zones** (the
  editor below; the catalog's zone properties, for ElastoPlasm c₀, cᵣ [kPa], φ [°], T, E [MPa], ν, ρ₀; φ hidden when the options pick
  `material.plastic = "VM"` (its `show_if`); "+ New zone ▾" opens a menu of shapes, `new_zone!(e, shape)`: Draw (a polygon,
  vertex by vertex), Rectangle (two opposite corners; every corner lying on another zone's edge splits
  it), Triangle (three corners, closes by itself), Circle (centre, then a point on it); the vertex last
  clicked or dragged (`picked`, outlined on the canvas) has one row, x and z typed exactly and "×"
  deleting it, `set_picked!`/`remove_picked!`, which move or delete it in every zone sharing it, through
  undo; the vertices are not listed), **Options** (one row per catalog option, pages by the first part of its key:
  checkbox, name, drop-down; numeric fields take any typed value plus presets; picking or typing a
  value ticks the box, unticking restores the default; multi-selections such as `plot.what` are
  checkbox grids; for ElastoPlasm every `get_option()` field but `dtype.T0`/`bits` (the one `precision` choice sets
  them) and `solution` (it follows the workflow)) and **Create** (named Run until the user pointed out that a generic
  layout has no run: the catalog's run settings under a Run heading when it has some, a drop-down or a field
  each, shown by `show_if`: for ElastoPlasm the workflow, then `te`/`tg`/`tep` or the load steps; the host's checks
  live, `check_case` for ElastoPlasm; then **Create…** and **Cancel**).
- A picture to draw on (a scanned section or a map, any image Qt reads), drawn under the grid at an
  adjustable opacity (`Backdrop`: file, pixel size reported by a hidden QML `Image`, lower-left corner and
  width in metres, the height from the aspect). The Domain step starts with a choice, *Domain* (sized by
  hand) or *Picture to draw on* (`use_picture!`, `editor.picture`): the latter opens the file dialog and
  shows the picture's own framed panel, in the order picture → scale → placement; *Domain* again hides the
  picture and keeps it. (It used to be a sub-section in the middle of the domain fields, with a separate
  "Apply" and "Fit domain to picture"; the user found the import hard to follow.) The picture opens across
  the domain's width; its scale is a scale bar ("Measure a scale bar": click its two ends on the canvas, a
  crosshair meanwhile, then type its length and Enter: `apply_scale!(e, text)`, a uniform scale about the
  picture's corner) or its width typed. **The domain follows the picture**: once it is loaded, and after
  every change of scale or lower-left corner (`x0`, `z0`, 0 by default; negative crops it), `fit_domain!`
  sets the domain from the origin to the picture's far corner, in whole square elements of the current `h`
  (the nearest number: a measured scale gives an arbitrary size, otherwise slightly elongated cells, which
  `check_case` flags under dynamic relaxation beyond 1%); `Lx`/`Lz` typed afterwards still change it. Checked
  offscreen on a 1000 × 500 px section with a 5 m scale bar (49 px/m): measured on the picture shown 10 m
  wide, it gives a picture 20.41 m wide and a 20.5 × 10 m domain of 41 × 20 elements at `h = 0.5`. A drawing
  aid only: not saved with the case (TOML).
  - **Loading, two Qt traps** (both hit by the user on a macOS screenshot, invisible in the tests until the
    window test loaded a picture, `EasyDraw.window(editor; shots)`): the hidden `Image` loads a local file at once,
    inside the refresh that set its source, so its size report (another refresh) is deferred with
    `Qt.callLater` (called directly: a binding loop on `source`; `asynchronous: true` delayed the report
    past the first paints instead). And the canvas's `loadImage` of a picture Qt already holds (the hidden
    `Image` loaded it) completes at once *without* `imageLoaded`, so the canvas draws it in the same pass
    when `isImageLoaded` turns true; before, the picture stayed blank until the next repaint (another tab, a zoom).
- Only the tabs and the status line are always visible. "Start from a file…" (top of the Domain
  step) loads a file through the host's `load` (`load!`; a case TOML for ElastoPlasm, a layout TOML by default;
  replaces the whole state and clears undo). "Load catalog…" beside it replaces the catalog with a TOML file's (asked
  for: a friend loads their own solver's catalog; `load_catalog!`): each zone keeps all its properties, shown or carried
  by the new catalog's list (in the undo steps too), still-valid option overrides and run settings kept, the old
  catalog's others dropped, a host's own run settings kept, an unknown side kind set to the default; not undoable. `draw`'s
  `catalog` may itself be the file (`Catalog(file)`). Checked by loading ElastoPlasm's catalog file
  (`write_catalog(file, get_catalog())`) into a window opened with a two-field catalog: the Options tab appears, a zone's
  `rho` stays carried. **Create…** (end of the
  Create step) runs `check_case` first (`ready!`: errors block it, warnings need a second click), then opens
  Qt's folder dialog (choose, or make with its New Folder button): `create!` sets the run's folder
  (`run.folder`), writes `<name>_case.toml` there and closes with the case;
  `setup_simulation` then writes the simulation (`.jld2`, plots, outputs) into that folder instead of
  `dump/<name>/` (a plain name, as for cases made in code). There is no Save…: the case file is written by
  Create…. The dialogs start in the current directory (`draw(…; startdir=pwd())`), not in `dump/`,
  the solver's scratch folder (`clean_dump()` empties it).
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
  one (`delete_vertex!`); each zone is drawn where it applies, since a point belongs to the last zone containing it:
  its fill and solid outline outside the later zones that may cover it (`later`, by bounding box; canvas clips,
  even–odd as `Polygon`'s rule), faint dashes under them (asked for: a first region crossed by later ones still showed
  its whole polygon, though its points already went to the later zones). **Drawing a zone cuts it out of the earlier
  zones** it overlaps (`cut_earlier!`, asked for next; the user chose automatic over on demand): their outlines then
  follow it, sharing its vertices (where outlines cross, the new zone gets a vertex too), so a zone may get several
  outlines (parts, holes); one covered entirely is removed; a cut circle becomes a polygon. The points do not change
  (the new zone won them anyway), checked against the priority rule in EasyDraw's tests, but for a cut circle's next to
  its arc; a circle cuts as an inscribed 32-gon, whose slivers it still wins. Deleting or moving the later zone then
  leaves a gap (undo restores it, one step with the drawing). `cut` (`cut.jl`): the edges of both split where they
  cross or touch, each piece kept when the difference lies on exactly one side of it (probed 1e-7·L away with the
  even–odd rule), chained into outlines turning most to the left; probing sides keeps snapped, shared edges exact.
  Overlaps remain possible (cases written in code, vertices dragged into another zone): the display above shows them.
  The selected zone stands out (asked for:
  a 3 px outline against 2 px was hard to tell): drawn last, a glow and a thick line inside its own area, also along
  the later zones' edges that cut into it, larger white-edged vertices and a stronger fill, the others thin and
  faded, its colour as a swatch beside the zone list; a click inside a zone selects
  it (the topmost, as `case_points` assigns; `clicked`, on release, unless the press became a pan); navigating (asked for once zoomed in to draw: there was no pan, and no way back but
  zooming out): a mouse wheel zooms about the pointer, a trackpad pans with two fingers (`pixelDelta`) and zooms
  with a pinch (`PinchHandler`) or Cmd/Ctrl + scroll (`wheel!`); dragging the empty canvas (a press that starts no
  edit, any tab) or the middle button (anywhere, also while drawing) pans, the grabbed point kept under the cursor
  (`grab`, `pan!`); "Fit view" (canvas corner) recentres on the domain (`fit!`, as switching tabs does; the same fit on every
  step: the Domain step used to show 1.4 times the domain, room to drag its corner, which left the canvas half
  blank; a drag continues outside the canvas and the view refits on release). x and z
  keep equal scales and the axes fill the canvas: `transform` widens the visible area about its centre along one
  axis, so there are no blank margins outside the axes, only coordinate space beside a domain whose proportions
  differ from the canvas (the user chose this over independent scales, which fill the canvas with the domain but
  stretch the picture, circles and slopes, and would need snapping in screen distances); undo/redo (Cmd/Ctrl+Z,
  Shift+Cmd/Ctrl+Z, buttons) keeps snapshots of the zones. The points to be created are previewed
  (`drawn_points`).

### Recording a run (HDF5 output)

```julia
case = slump_case(L, nel; cli()..., plot=(; status=false),
                  output=(; status=true, every=0.5, what=(; point=(; solid=["u", "p", "epII"]))))
jld2 = setup_simulation(case)
execute!(jld2; workflows=[elastodynamic!, elastoplastic!])         # writes <name>_output.h5 beside the .jld2
f    = replace(jld2, "_simulation.jld2" => "_output.h5")
output_series(f, "point/solid/p"; point=30)                         # (; t, values, loadstep)
plot_output(f, "point/solid/p"; frame=7); animate_output(f, "point/solid/p"); plot_series(f, "point/solid/p"; point=30)
```

- **Configuration:** the `output` solver section, `(; status, every, what)`; `what` mirrors
  the data model, `point = (; base, solid, fluid, thermal)` and `mesh = (; solid, thermal)`, each a
  list of keys from `quantities()` (`src/home/api/io/quantities.jl`, shaped like `what`: point base
  `x J n Omega ell`, solid `u v p dev epII epV rho q umag vmag c0 phi0`, thermal `T`; mesh solid `v a m`,
  thermal `T`). `quantities()` is the one table of what can be plotted or recorded: per quantity its
  extractor (a frame in SI), readable `name`, SI `unit`, and for the plottable point ones (those with a
  `label`) the displayed unit and factor, colour map and limits. `plot.what` picks from the plottable
  keys (`plottable()`), `get_option()` derives both lists from it, the window and the CLI show the
  readable names (`shown_name`; an older case's `plot.what` written by name still loads). A new quantity
  is one entry there. Partial selections merge per field; a phase the problem lacks is skipped with a
  warning. Multi-selections are any non-numeric option whose default is a list (`multi_select`):
  `plot.what` and these lists, checkbox grids in the window, all on its "output" page (`point · solid`, `mesh · solid`, …; a list
  with nothing to pick, such as the fluid phase's, is not shown).
- **File:** one array per field, frames last (`/point/<phase>/<field>`, `/mesh/<phase>/<field>`),
  chunked per frame, deflate 3, Float64; `/frames/t` (a time, or a load step:
  `/frames/loadstep` 1), `/frames/step`, `/frames/workflow` (the workflow's name); the grid,
  the initial point properties and attributes (case TOML, version, laws). Each field's dataset has the
  attributes `name`, `unit` (SI) and, when plottable, `label` (LaTeX); `/frames/t` has `unit = "s"`
  (`load_output(f).info`, `output_frame`'s `info`). Plain HDF5 (h5py etc.).
- **When:** explicit workflows stop on the `plot.freq` grid and the `output.every` grid, each ending at
  the stage's end time (`checkpoints`, `explicit/get.jl`; without output exactly the plot stops, so
  results are bit-identical with output off); `elastoquasistatic!` records each load step; `execute!`
  records the initial state and the state each workflow ends on. The recorder lives in
  `task_local_storage` for the run (`record!` does nothing without one), so workflow signatures and
  closures are unchanged. The file is `output` (default `<folder>/<name>_output.h5`, `output_file(sim)`):
  `execute!(jld2)` appends (it continues the saved state), `execute` and an in-memory
  `execute!(sim; output, append=false)` start a new one, `setup_simulation` removes a stale one. Adding stops changes the time-step sequence
  (`get_dt` clips at each stop), so a run with a finer `every` than `plot.freq` differs slightly
  from one without output.
- **After the run:** `output_series`, `output_frame`, `load_output` read it; `plot_output`,
  `animate_output`, `plot_series` draw from the file alone (live-plot colour maps, units and limits
  where the field matches a plot variable). Checked: a frame drawn from the file matches the live
  plot of the same time; `test_output.jl` covers shapes, cadence, last frame = saved state, append,
  skipped phases, load steps and the plotting.

### Running on a GPU (explicit path)

```julia
using CUDA, ElastoPlasm               # CUDA in the environment: UnifiedBackend's CUDA extension registers the GPU
case = slump_case(L, nel; backend=(; select="device"))
execute!(setup_simulation(case))      # elastodynamic!/elastoplastic! on the GPU; results back on the host
```

- `backend.select = "device"` picks the first GPU UnifiedBackend registered (CUDA, ROCm); without one, the host, with
  a message. The kernels are built on its backend (`exec_backend`), the workflows move the problem there once and
  copy it back at the checkpoints that plot or record, and at the end (`to_device`/`to_host!`, `core/device.jl`), so
  plots, HDF5 output and the `.jld2` are as on the host; a saved simulation reloads anywhere, on that machine's device
  (`persistence.md`). Same precision as on the host (`dtype`); Float64 atomics need compute capability ≥ 6.0, and
  consumer GPUs run Float64 much slower than Float32.
- Not yet on a device: dynamic relaxation (`elastoquasistatic!`) and thermal (`thermodynamic!`) refuse it with an
  `ArgumentError`.
- Results differ from a host run by round-off (`gotchas.md`): compare with a relative tolerance.
- **Validation on an NVIDIA machine:** `julia -t 1 test/cuda/validate_cuda.jl` from the repository root (cloned with
  `--recursive`) sets up its own environment (`test/cuda`, gitignored), reports the GPU, compiles and runs a few steps
  of every configuration (all errors collected), compares device and host runs, times a step on both and checks that a
  file saved there loads without CUDA; it writes `dump/cuda_validation/report.toml`. `EP_DEVICE=jlarrays` runs it with
  JLArrays instead, without a GPU (small sizes).

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
trip, `check_case`/`estimate_dt` (= `get_dt`, 2D and 3D), a case read from its layout and settings (a copy, 3D sides, one geometry with two settings files), zone constants reaching the particles,
free option values, workflows, a short `execute!`, and the `Simulation` (build, save, load, in memory equal to the file run, older files).
`test_easydraw.jl` runs EasyDraw's own tests (`lib/EasyDraw/test/runtests.jl`, also `] test EasyDraw` in an
environment that has it): shapes, zones, layouts and catalogs as TOML, sampling (`zone_at`,
`property`, `points`), the partition (areas by hand, shared edges, sides, holes, a circle) and its `.geo`, the editor
with a small test catalog and with none, driven as a user would (drawing, snapping and shared edges, drag, undo,
cutting, circles, the picture, domain, sides, zone properties, options, run, Create… (its warning confirmation, the
folder) and loading a file, the caches and the readout, the JSON the view reads); then the window offscreen, a PNG
of each step, twice in one session, then once on a picture (a BMP written by
the test: Qt measures it, the domain follows). `test_draw.jl` covers ElastoPlasm's side: `get_catalog` (every
solver option, the zone constants and `show_if`, the run fields, as TOML), a case's layout through the window and back unchanged, a collision's
preview points equal to the problem's, a case drawn through the editor (options, run, `check_case`, Create… writing a
loadable `<name>_case.toml`, loading a case file), and `draw_case`'s window offscreen.

`test_device.jl` runs the explicit workflows on a device: JLArrays (GPUArrays' reference array, on the CPU; a
test-only dependency, `[extras]`, installed by `Pkg.test()`; the file is skipped without it), registered as a
UnifiedBackend device, with host element access forbidden (`allowscalar(false)`). Five 2D configurations, 3D and
Float32 match the host run to round-off (Float64 1e-9 relative, Float32 1e-5; integers and Bools exactly; see
`gotchas.md`), the device copy holds device arrays only and returns to the same types, and dynamic relaxation and
thermal refuse a device. To run it outside `Pkg.test()`, put an environment with JLArrays on the load path
(`JULIA_LOAD_PATH="@:<env with JLArrays>:@stdlib"`).

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
config merged by `solver_instr` — this is the single place that defines every tunable
knob and its out-of-the-box value. Its physics sections (`material`, `grf`, `bcs`) are a case's (`case.laws`,
`case.grf`, `case.bcs`): `Solver` refuses them, and the solver holds the others.

- `solution` — `"explicit"`/`"implicit"`, the solver's `solution` marker (`ExplicitSolution()`/`ImplicitSolution()`,
  `get_solution`, `SOLUTIONS`): one `Solver` struct, `ExplicitSolver`/`ImplicitSolver` its aliases. Printed (as
  `solution_name(solver)`) in `execute.jl`'s generated filenames and `logs.jl`'s startup banner. Picking `solution=
  "explicit"` then calling an implicit-only workflow fails with `MethodError` at the
  workflow call, not silently — `elastodynamic!`/`elastoplastic!` take an
  `ExplicitSolver` specifically. `ignite()` is the one place taking any `Solver`
  instead, since both paths call it regardless of `solution`.
- `dtype` — arithmetic precision (`bits`, element types `T0`)
- `basis` — `which` (`"bsmpm"`/`"gimpm"`/`"smpm"`/`"mlsmpm"`, see `get_basis`), `trsfr` (P2G/G2P transfer scheme: `"std"`/`"tpic"`/
  `"apic"`, see `get_transfer`), `C_pf` (PIC/FLIP blend), `domain` (GIMP's particle-domain update, read by
  `gimpm` only: `"stretch"`, the default, `ℓ = ℓ₀ ∘ diag(U)` after Charlton et al. 2017 and AMPLE Eq 34, or
  `"undeformed"`, uGIMP, `ℓ = ℓ₀`; see `get_domain`, `particle_domain`). `trsfr`/`C_pf` live here
  rather than a separate `transfer` section since `Basis` (the struct) owns both
  `kind` and `transfer` as sibling fields — see "Transfer scheme dispatch" in
  `architecture.md`.
- `stab` — `locking` (F-bar volumetric locking correction on/off), `damping`, `musl`
  (MUSL velocity reprojection on/off — lives here rather than under `basis`/`transfer`
  since it's a stabilization technique applied regardless of transfer scheme).
- `bcs` — `dirichlet` boundary condition matrix, one `[lower upper]` row per dimension (a case's: its sides; given
  among a case constructor's solver options, it moves to the case)
- `grf` — Gaussian random field for heterogeneous initial cohesion (`status` toggles it on;
  `covariance` `"gaussian"` or `"exponential"`; `param`: standard deviation `σ`, per-axis correlation lengths
  `Iₓ`, harmonics `Nₕ`, max wave number `kₘ`). Evaluated at the material points' own coordinates
  by `property_field`/`random_field` (`init/mpts/GRF.jl`), following Räss et al. (2019), in
  `refs/ElastoPlasm/grf/`.
- `material` — `elastic` (`"hencky"`/`"improved_hencky"`/`"hypoelastic"` — picks both
  the strain kinematics and the elastic law, see `AbstractElasticLaw` in `lagrangian.jl`;
  `"hypoelastic"` is the small-strain Jaumann-rate path, the two Hencky laws are finite
  strain) and `plastic` (`"DP"`/`"VM"` — the only two `setup_cmp` builds; any other value
  is an `ArgumentError`; `ELASTIC_LAWS`/`PLASTIC_LAWS` map the names to the points' types; the `retmap` kernel
  dispatch lives on `Point`'s `CM` type parameter, not a runtime string — see "DP/J2
  retmap kernel unification" in `planned-improvements.md`)
- `nonloc` — non-local plastic strain regularization (`status`, `ls` length scale)
- `plot` — `status`, `freq` (plot every N `Time` checkpoints), `dpi`, `what` (keys of the plottable
  quantities, `plottable()` in `quantities.jl`, e.g. `["p", "epII"]`)
- `backend` — `select` (`"host"`, or `"device"`: a GPU UnifiedBackend registered, through its CUDA/ROCm extension,
  else the host with a message); `build_solver` stores the selected UnifiedBackend `Platform` as `backend.exec`. (A
  `backend.distributed` option existed for multi-device selection, which nothing used; removed with UnifiedBackend 0.1.0,
  so a case file setting it is refused as an unknown option.)

Two ways to change solver behaviour:

1. **Per-run override (preferred, no code change)** — pass any of these keys as
   keyword arguments to a case constructor (`slump_case`, every section) or to `Solver(case; …)`/`Solver(; dim, …)`
   (the numerical sections); `solver_instr` only keeps
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
struct itself (like `run` did), it also needs threading through
`Solver`'s field list (`src/boot/needs/types/solver.jl`), its partial constructor, `solver_config`/`stored_solver`
and the corresponding positional arg in `api/solver/build.jl`'s final constructor call. A new physics option (one a case
holds, like the laws) goes into `PHYSICS_SECTIONS` and a `Case` field instead.

### Using `cli()`

`cli(; ui::Bool=false)` turns `get_default()`'s config tree into a `Dict{Any,Any}` of
kwargs suitable for splatting into `slump_case` as `cli()...` (not into `Solver(; …)`, which refuses its physics
sections).

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
  `get_option().basis.trsfr` → `["std", "tpic", "apic"]`). `backend.select = "device"` runs the
  explicit workflows on a GPU, while dynamic relaxation and thermal refuse it (see "Running on a
  GPU" above). It offered `"cuda"`/`"rocm"` until UnifiedBackend's adoption, values the selection
  refused (`ArgumentError`).
