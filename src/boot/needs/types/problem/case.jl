# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Simulation case: the declarative description a problem is built from (domain, zones, materials,
# boundaries, solver overrides, run settings), built into a `Problem` (`Problem(case)`, `setup_simulation`)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractShape, Polygon, Block, Ball, Region, Zone, Case, PointData, MATERIAL, inpolygon

"""
    MATERIAL

The default material of a `Case`, in SI units with angles in radians: Young's modulus `E`
[Pa], Poisson's ratio `nu`, grain density `rho0` [kg/m³], cohesion `c0` and residual cohesion `cr`
[Pa], friction angle `phi` [rad], softening modulus `Hp` [Pa], initial porosity `n0`, temperature
`T` [K], specific heat capacity `cp` [J/(kg·K)], thermal conductivity `k` [W/(m·K)].
"""
const MATERIAL = (; E=1.0e6, nu=0.3, rho0=2700.0, c0=20.0e3, cr=4.0e3, phi=20.0*π/180, Hp=-30.0e3, n0=0.1,
                    T=293.15, cp=10.0, k=3000.0)
# what a zone may set: any material constant, and an initial velocity `v` [m/s]
const ZONE_PROPS = (keys(MATERIAL)..., :v)
# the run settings of a case: output folder under dump/, times [s] for `setup_time` (elastic stage,
# gravity ramp, plastic stage), workflow (`WORKFLOWS`), load steps of the quasi-static workflows
const RUN = (; folder="problem", te=10.0, tg=10.0, tep=5.0, workflow="elastodynamic! + elastoplastic!", lstps=10)
const BOUNDARY_KINDS = (:free, :roller, :fixed)
# the shapes (`Polygon`, `Block`, `Ball`, `Region`, their `inside` rule) and `Zone` are EasyDraw's (lib/EasyDraw),
# imported in boot.jl; a case's zones may set only `ZONE_PROPS`, checked by `Case`

"""
    Case

A simulation case: a model as an EasyDraw `Layout`, read as ElastoPlasm needs it. The layout is the data (domain,
zones, a boundary kind per side, the material's defaults, the solver options, the run settings), drawn in the window
or written in code; the case checks it and holds its typed reading, which `Problem(case)` and `setup_simulation`
build a problem from. A case keeps its own copy of the layout: to change a model, change a layout and build a new case.

# Fields
- `layout::EasyDraw.Layout`: the model: `L` [m] and `nel` (2D x–z, or 3D x, y, z; the lower corner at the origin),
  `ni`×`ni` material points per element, the `zones` (`Zone`, a point belongs to the last zone containing it), and
  `sides`, `defaults`, `options`, `run`, read below. Its file is the case's file (`save_case`, `load_case`).
- `material::NamedTuple`: the material constants: `MATERIAL`, with the layout's `defaults`.
- `bcs::Matrix{Symbol}`: the boundary condition of each side, one `[lower upper]` row per axis (x, then y in 3D,
  then z): the layout's `sides`, `"left"`/`"right"`, `"front"`/`"back"`, `"bottom"`/`"top"`, each `:free`, `:roller`
  (a side not set) or `:fixed`.
- `solver::NamedTuple`: the solver overrides: the layout's `options` as values (`option_values`), merged per field
  over `get_default()` by `get_solver`; `solution` follows the workflow unless set.
- `thermal::Bool`: whether to build the thermal phase: the layout's run setting `thermal`, false when absent.
- `run::NamedTuple`: the run settings, `RUN` with the layout's (`case_run`: `folder`, `te`, `tg`, `tep`, `workflow`,
  `lstps`). The folder is the layout's `run["folder"]`, else the layout's own folder (a file's), else `RUN.folder`.

# Constructors
    Case(layout)
    Case(; L, nel, ni=2, zones=[], material=(;), bcs=nothing, solver=(;), thermal=false, run=(;))

The keyword form builds the layout, then reads it: `bcs` (or `solver.bcs.dirichlet`, as the solver options give them;
a row per axis, extra rows ignored) as the sides, `material` (overrides of `MATERIAL`) as the defaults, `solver`
(overrides of `get_default()`) as the options in their file form (`option_labels`), and every run setting, with
`thermal`.

# Throws
- `ArgumentError`: when a zone sets something other than a material constant or a velocity `v`, a default is not a
  material constant, a side's kind is not one of `BOUNDARY_KINDS`, an option or a run setting is unknown or invalid,
  or the layout itself is invalid.
"""
struct Case
    layout   :: EasyDraw.Layout
    material :: NamedTuple
    bcs      :: Matrix{Symbol}
    solver   :: NamedTuple
    thermal  :: Bool
    run      :: NamedTuple
end
function Case(layout::EasyDraw.Layout)
    l = deepcopy(layout)                                     # the case's own: nothing outside can change it
    check_zone_props(l.zones)
    run = Dict{String,Any}(l.run)
    thermal = Bool(pop!(run, "thermal", false))
    haskey(run, "folder") || isempty(l.folder) || (run["folder"] = l.folder)
    return Case(l, case_material(l.defaults), case_bcs(l.sides, length(l.L)), case_solver(l.options), thermal, case_run(run))
end
function Case(; L, nel, ni::Integer=2, zones=Zone[], material=(;), bcs=nothing, solver=(;), thermal::Bool=false, run=(;))
    D = length(L)
    s = named_tuple(solver)
    if haskey(s, :bcs)                                       # boundary conditions passed with the solver options
        isnothing(bcs) && (bcs = s.bcs.dirichlet)
        s = Base.structdiff(s, NamedTuple{(:bcs,)})
    end
    b = isnothing(bcs) ? fill(:roller, D, 2) : Symbol.(bcs)
    size(b, 1) > D && size(b, 2) == 2 && (b = b[1:D, :])    # more rows than axes (the solver default has three)
    size(b) == (D, 2) || throw(ArgumentError("bcs is one [lower upper] row per axis ($(D)×2)"))
    sides    = Dict(side => string(b[axis, k]) for (side, (axis, k)) ∈ side_axes(D))
    defaults = Dict{String,Any}(string(k) => v for (k, v) ∈ pairs(named_tuple(material)))
    settings = config_dict(case_run(run))
    thermal && (settings["thermal"] = true)
    layout = EasyDraw.Layout(; L, nel, ni, zones, sides, defaults, options=option_labels(config_dict(s)), run=settings)
    return Case(layout)
end

# the sides of a D-dimensional domain, each as its axis and end (1: lower, 2: upper) in `bcs`
function side_axes(D::Integer)
    if D == 2
        return ("left" => (1, 1), "right" => (1, 2), "bottom" => (2, 1), "top" => (2, 2))
    end
    return ("left" => (1, 1), "right" => (1, 2), "front" => (2, 1), "back" => (2, 2), "bottom" => (3, 1), "top" => (3, 2))
end
# a Dict of sections (from a file or the window) as a NamedTuple, and a NamedTuple as it is
function named_tuple(x)
    return x isa AbstractDict ? (; (Symbol(k) => named_tuple(v) for (k, v) ∈ x)...) : x
end
# what a zone may set: material constants (numbers) and a velocity (numbers)
function check_zone_props(zones)
    for z ∈ zones
        bad = setdiff(keys(z.props), ZONE_PROPS)
        isempty(bad) || throw(ArgumentError("zone \"$(z.name)\": unknown properties $(collect(bad)); known: $(ZONE_PROPS)"))
        all(((k, v),) -> k == :v ? v isa Vector{Float64} : v isa Float64, z.props) ||
            throw(ArgumentError("zone \"$(z.name)\": its material constants are numbers, its velocity v a list of numbers"))
    end
    return nothing
end
# the material: MATERIAL with the layout's defaults
function case_material(defaults)
    m   = Dict{Symbol,Any}(Symbol(k) => v for (k, v) ∈ defaults)
    bad = setdiff(keys(m), keys(MATERIAL))
    isempty(bad) || throw(ArgumentError("unknown material constants $(collect(bad)); known: $(keys(MATERIAL))"))
    return (; (k => Float64(get(m, k, getfield(MATERIAL, k))) for k ∈ keys(MATERIAL))...)
end
# the boundary conditions: the layout's sides, roller where a side is not set
function case_bcs(sides, D)
    bcs = fill(:roller, D, 2)
    for (side, (axis, k)) ∈ side_axes(D)
        kind = Symbol(get(sides, side, "roller"))
        kind ∈ BOUNDARY_KINDS || throw(ArgumentError("side $side: $kind is not one of $(BOUNDARY_KINDS)"))
        bcs[axis, k] = kind
    end
    return bcs
end
# the solver overrides: the layout's options as values, checked against get_default() (an unknown field errors)
function case_solver(options)
    s = named_tuple(option_values(options))
    merge_config(get_default(), s)
    return s
end

"""
    case_run(run) -> NamedTuple

Run settings (`RUN`: `folder`, where the simulation is written, a name under `dump/` or a path;
`te`, `tg`, `tep`, `workflow`, `lstps`), given as a `Dict` or `NamedTuple`, missing ones from the
defaults; `workflow` is one of `WORKFLOWS`. The `fid`/`dir` of older cases are read as `folder`.
"""
function case_run(run)
    unknown = setdiff(string.(keys(run)), [string.(keys(RUN))..., "fid", "dir"])
    isempty(unknown) || throw(ArgumentError("unknown run settings $(unknown); known: $(string.(keys(RUN)))"))
    val(k, d=getfield(RUN, k)) = run isa AbstractDict ? get(run, string(k), get(run, k, d)) : get(run, k, d)
    old = (d = val(:dir, ""); isempty(d) ? val(:fid, RUN.folder) : d)  # older cases
    r = (; folder=String(val(:folder, old)), te=Float64(val(:te)), tg=Float64(val(:tg)), tep=Float64(val(:tep)),
           workflow=String(val(:workflow)), lstps=Int(val(:lstps)))
    r.workflow ∈ WORKFLOWS || throw(ArgumentError("run.workflow = \"$(r.workflow)\" is not one of $(WORKFLOWS)"))
    (r.te ≥ 0 && r.tg ≥ 0 && r.tep ≥ 0 && r.lstps ≥ 1) || throw(ArgumentError("run times must be ≥ 0 and lstps ≥ 1"))
    return r
end
# the workflows a case may run (`case_workflows`); the quasi-static ones need `solution = "implicit"`
const WORKFLOWS = ("elastodynamic! + elastoplastic!", "elastodynamic!", "elastoplastic!",
                   "elastoquasistatic!", "elastoquasistatic! (plastic)")
quasistatic(workflow) = startswith(workflow, "elastoquasistatic!")


"""
    PointData

The material points a problem starts from, one column or entry per point: positions `xp` (D×nmp)
and initial velocities `vp`, points per element `ni`, the `zone` each belongs to (0 for uniform points),
and their constants (`c0`, `cr` [Pa], `phi` [rad], `Hp` [Pa], `E` [Pa],
`nu`, `rho0` [kg/m³], `n0`, `T` [K], `cp`, `k`). Built by `case_points` (or uniform,
`PointData(xp, ni)`), read by `setup_mpts`.
"""
struct PointData
    xp   :: Matrix{Float64}
    vp   :: Matrix{Float64}
    ni   :: Int
    zone :: Vector{Int}
    c0   :: Vector{Float64}
    cr   :: Vector{Float64}
    phi  :: Vector{Float64}
    Hp   :: Vector{Float64}
    E    :: Vector{Float64}
    nu   :: Vector{Float64}
    rho0 :: Vector{Float64}
    n0   :: Vector{Float64}
    T    :: Vector{Float64}
    cp   :: Vector{Float64}
    k    :: Vector{Float64}
end

"""
    PointData(xp, ni; material=MATERIAL) -> PointData

Points at `xp` (D×nmp, `ni` per element and axis), at rest, all of one `material`.
"""
function PointData(xp::AbstractMatrix, ni::Integer; material=MATERIAL)
    n = size(xp, 2)
    f(k) = fill(Float64(getfield(material, k)), n)
    return PointData(Float64.(xp), zeros(size(xp)), Int(ni), zeros(Int, n), f(:c0), f(:cr), f(:phi), f(:Hp), f(:E), f(:nu),
                     f(:rho0), f(:n0), f(:T), f(:cp), f(:k))
end
