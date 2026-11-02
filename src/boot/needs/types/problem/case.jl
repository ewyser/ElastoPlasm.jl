# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Simulation case: the declarative description a problem is built from (domain, zones, materials,
# boundaries, solver overrides, run settings), built with a `Solver` into a `Problem` (`Problem(case, solver)`)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractShape, Polygon, Block, Ball, Region, Zone, Layout, Settings, load_layout, load_settings, save_settings,
       Case, PointData, get_material, inpolygon

"""
    get_material() -> NamedTuple

The default material of a `Case`, in SI units with angles in radians: Young's modulus `E` [Pa], Poisson's ratio `nu`,
grain density `rho0` [kg/m³], cohesion `c0` and residual cohesion `cr` [Pa], friction angle `phi` [rad], softening
modulus `Hp` [Pa], initial porosity `n0`, temperature `T` [K], specific heat capacity `cp` [J/(kg·K)], thermal
conductivity `k` [W/(m·K)]. A case's settings override it (`defaults`, then each zone's values).
"""
function get_material()
    return (; E=1.0e6, nu=0.3, rho0=2700.0, c0=20.0e3, cr=4.0e3, phi=20.0*π/180, Hp=-30.0e3, n0=0.1, T=293.15,
              cp=10.0, k=3000.0)
end
# what a zone may set: any material constant, and an initial velocity `v` [m/s]
function get_zone_properties()
    return (keys(get_material())..., :v)
end
# the default run settings of a case: output folder under dump/, times [s] for `setup_time` (elastic stage, gravity
# ramp, plastic stage), workflow (`get_workflows()`), load steps of the quasi-static workflows
function get_run()
    return (; folder="problem", te=10.0, tg=10.0, tep=5.0, workflow="elastodynamic! + elastoplastic!", lstps=10)
end
# the boundary condition a side may have
function get_boundary_kinds()
    return (:free, :roller, :fixed)
end
# the shapes (`Polygon`, `Block`, `Ball`, `Region`, their `inside` rule) and `Zone` are EasyDraw's (lib/EasyDraw),
# imported in boot.jl; a case's zones may set only `get_zone_properties()`, checked by `Case`

"""
    Case

A simulation case: a geometry (an EasyDraw `Layout`) with its physical properties assigned, from its settings (an
EasyDraw `Settings`): what is solved. The layout is where things are (domain, grid, zones by name); the settings hold
what they are (each zone's material, the material's defaults, a boundary kind per side, the constitutive laws) and how
a solver built on the case solves it (the numerical options, the run settings). The case checks all of them and holds
the typed reading of its physics; `Solver(case)` reads the numerics, `Problem(case, solver)` builds what both
describe. A case keeps its own copy of both: to change a model, change a layout or settings and build a new case. One
geometry can take several settings, `Case(layout, settings_a)` and `Case(layout, settings_b)`.

# Fields
- `layout::EasyDraw.Layout`: the geometry: `L` [m] and `nel` (2D x–z, or 3D x, y, z; the lower corner at the
  origin), `ni`×`ni` material points per element, the `zones` (a name and a shape each; a point belongs to the last
  zone containing it).
- `settings::EasyDraw.Settings`: the values, as plain data: `zones` (by zone name), `defaults`, `sides`, `options`,
  `run`. The case reads the physics below and its folder (`case_folder`); `Solver(case)` the numerical options and
  the run settings. With the layout, the case's file (`save_case`, `load_case`).
- `material::NamedTuple`: the material constants: `get_material()`, with the settings' `defaults`.
- `properties::Dict{String,NamedTuple}`: each zone's own constants, by zone name: `get_material()` keys, or an initial
  velocity `v`; where a zone does not set one, the material's.
- `bcs::Matrix{Symbol}`: the boundary condition of each side, one `[lower upper]` row per axis (x, then y in 3D,
  then z): the settings' `sides`, `"left"`/`"right"`, `"front"`/`"back"`, `"bottom"`/`"top"`, each `:free`,
  `:roller` (a side not set) or `:fixed`.
- `laws::NamedTuple`: the constitutive laws, `(; elastic, plastic)`: the settings' options `material.elastic`
  (`"hencky"`, `"improved_hencky"` or `"hypoelastic"`) and `material.plastic` (`"DP"` or `"VM"`), `get_default()`'s
  where not set.
- `grf::NamedTuple`: the random field of the cohesion, `(; status, covariance, param)`: the settings' options `grf`,
  `get_default()`'s where not set.
- `thermal::Bool`: whether to build the thermal phase: the run setting `thermal`, false when absent.

# Constructors
    Case(layout, settings)
    Case(; L, nel, ni=2, zones=[], properties=Dict(), material=(;), bcs=nothing, solver=(;), thermal=false, run=(;))

The keyword form builds both, then reads them: `zones` (a name and a shape each) as the layout; `properties` (by zone
name, e.g. `Dict("rock" => (; E=5e7))`) as the zones' values; `material` (overrides of `get_material()`) as the
defaults; `bcs` (or `solver.bcs.dirichlet`, as the solver options give them; a row per axis, extra rows ignored) as the
sides; `solver` (overrides of `get_default()`) as the options in their file form (`option_labels`); and every run
setting, with `thermal`.

# Throws
- `ArgumentError`: when a zone sets something other than a material constant or a velocity `v`, a default is not a
  material constant, a side's kind is not one of `get_boundary_kinds()`, an option or a run setting is unknown or
  invalid (all are checked here, also the numerical ones), or the layout and settings are invalid or do not agree
  (`EasyDraw.check_settings`).
"""
struct Case
    layout     ::EasyDraw.Layout
    settings   ::EasyDraw.Settings
    material   ::NamedTuple
    properties ::Dict{String,NamedTuple}
    bcs        ::Matrix{Symbol}
    laws       ::NamedTuple
    grf        ::NamedTuple
    thermal    ::Bool
end
function Case(layout::EasyDraw.Layout, settings::EasyDraw.Settings)
    l, s = deepcopy(layout), deepcopy(settings)              # the case's own: nothing outside can change them
    EasyDraw.check_settings(l, s)
    options = case_options(s.options)                        # every option, checked; the case keeps its physics
    case_run(l, s)                                           # every run setting, checked; the solver reads them
    d = get_default()
    laws = haskey(options, :material) ? merge_config(d.material, options.material; path="material.") : d.material
    grf  = haskey(options, :grf) ? merge_config(d.grf, options.grf; path="grf.") : d.grf
    return Case(l, s, case_material(s.defaults), case_properties(s.zones), case_bcs(s.sides, length(l.L)), laws, grf,
                Bool(get(s.run, "thermal", false)))
end
function Case(; L, nel, ni::Integer=2, zones=Zone[], properties=Dict(), material=(;), bcs=nothing, solver=(;),
                thermal::Bool=false, run=(;))
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
    layout = EasyDraw.Layout(; L, nel, ni, zones)
    return Case(layout, EasyDraw.Settings(; zones=properties, defaults, sides, options=option_labels(config_dict(s)), run=settings))
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
# each zone's constants, by zone name: material constants (numbers) and a velocity (numbers)
function case_properties(zones)
    out, known = Dict{String,NamedTuple}(), get_zone_properties()
    for (name, values) ∈ zones
        p   = Dict{Symbol,Any}(Symbol(k) => v for (k, v) ∈ values)
        bad = setdiff(keys(p), known)
        isempty(bad) || throw(ArgumentError("zone \"$name\": unknown properties $(collect(bad)); known: $(known)"))
        all(((k, v),) -> k == :v ? v isa Vector{Float64} : v isa Float64, p) ||
            throw(ArgumentError("zone \"$name\": its material constants are numbers, its velocity v a list of numbers"))
        out[name] = (; p...)
    end
    return out
end
# the material: `get_material()` with the layout's defaults
function case_material(defaults)
    m, d = Dict{Symbol,Any}(Symbol(k) => v for (k, v) ∈ defaults), get_material()
    bad  = setdiff(keys(m), keys(d))
    isempty(bad) || throw(ArgumentError("unknown material constants $(collect(bad)); known: $(keys(d))"))
    return (; (k => Float64(get(m, k, getfield(d, k))) for k ∈ keys(d))...)
end
# the boundary conditions: the layout's sides, roller where a side is not set
function case_bcs(sides, D)
    bcs, kinds = fill(:roller, D, 2), get_boundary_kinds()
    for (side, (axis, k)) ∈ side_axes(D)
        kind = Symbol(get(sides, side, "roller"))
        kind ∈ kinds || throw(ArgumentError("side $side: $kind is not one of $(kinds)"))
        bcs[axis, k] = kind
    end
    return bcs
end
# the option sections that are physics, read by the case; the others are numerics, read by `Solver(case)`
function get_physics_sections()
    return (:material, :grf)
end
# a case's options as values, checked against get_default() (an unknown field errors)
function case_options(options)
    s = named_tuple(option_values(options))
    merge_config(get_default(), s)
    return s
end

# the run settings of a layout's settings: `get_run()` with them (`thermal` aside); the folder, when they set none, the
# layout's own (a file's)
function case_run(l::EasyDraw.Layout, s::EasyDraw.Settings)
    run = Dict{String,Any}(s.run)
    delete!(run, "thermal")
    haskey(run, "folder") || isempty(l.folder) || (run["folder"] = l.folder)
    return case_run(run)
end
# where a case's simulation is written: its run folder, else its layout's file's folder, else `get_run().folder`
function case_folder(case::Case)
    return case_run(case.layout, case.settings).folder
end
# the run settings a solver holds: `case_run`'s without the folder, which is the case's (`case_folder`)
function solver_run(run)
    r = case_run(run)
    return (; te=r.te, tg=r.tg, tep=r.tep, workflow=r.workflow, lstps=r.lstps)
end

"""
    case_run(run) -> NamedTuple

Run settings (`get_run()`: `folder`, where the simulation is written, a name under `dump/` or a path;
`te`, `tg`, `tep`, `workflow`, `lstps`), given as a `Dict` or `NamedTuple`, missing ones from the
defaults; `workflow` is one of `get_workflows()`. The `fid`/`dir` of older cases are read as `folder`.
"""
function case_run(run)
    default, workflows = get_run(), get_workflows()
    unknown = setdiff(string.(keys(run)), [string.(keys(default))..., "fid", "dir"])
    isempty(unknown) || throw(ArgumentError("unknown run settings $(unknown); known: $(string.(keys(default)))"))
    val(k, d=getfield(default, k)) = run isa AbstractDict ? get(run, string(k), get(run, k, d)) : get(run, k, d)
    old = (d = val(:dir, ""); isempty(d) ? val(:fid, default.folder) : d)  # older cases
    r = (; folder=String(val(:folder, old)), te=Float64(val(:te)), tg=Float64(val(:tg)), tep=Float64(val(:tep)),
           workflow=String(val(:workflow)), lstps=Int(val(:lstps)))
    r.workflow ∈ workflows || throw(ArgumentError("run.workflow = \"$(r.workflow)\" is not one of $(workflows)"))
    (r.te ≥ 0 && r.tg ≥ 0 && r.tep ≥ 0 && r.lstps ≥ 1) || throw(ArgumentError("run times must be ≥ 0 and lstps ≥ 1"))
    return r
end
# the workflows a case may run (`case_workflows`); the quasi-static ones need `solution = "implicit"`
function get_workflows()
    return ("elastodynamic! + elastoplastic!", "elastodynamic!", "elastoplastic!", "elastoquasistatic!",
            "elastoquasistatic! (plastic)")
end
quasistatic(workflow) = startswith(workflow, "elastoquasistatic!")


"""
    PointData

The material points a problem starts from, one column or entry per point: the explicit contract between a case
(`case_points`) and the points built from it (`setup_mpts`, `build_solid_phase`, `build_thermal_phase`, `setup_cmp`).

# Fields
- `xp::Matrix{Float64}`: the positions, `D × nmp` [m].
- `vp::Matrix{Float64}`: the initial velocities, `D × nmp` [m/s].
- `ni::Int`: the points per element and axis.
- `zone::Vector{Int}`: the zone each point belongs to, its index in the layout (0 for uniform points).
- `c0`, `cr::Vector{Float64}`: the initial and residual cohesion [Pa].
- `phi::Vector{Float64}`: the friction angle [rad].
- `Hp::Vector{Float64}`: the softening modulus [Pa].
- `E::Vector{Float64}`: Young's modulus [Pa]; `nu`: Poisson's ratio.
- `rho0::Vector{Float64}`: the grain density [kg/m³]; `n0`: the initial porosity.
- `T::Vector{Float64}`: the temperature [K]; `cp`: the specific heat capacity [J/(kg·K)]; `k`: the thermal
  conductivity [W/(m·K)].

# Constructors
    PointData(xp, vp, ni, zone, c0, cr, phi, Hp, E, nu, rho0, n0, T, cp, k)
    PointData(xp, ni; material=get_material())

The second, points at rest of one material (`sfrac.jl`'s quadrature points).
"""
struct PointData
    xp   ::Matrix{Float64}
    vp   ::Matrix{Float64}
    ni   ::Int
    zone ::Vector{Int}
    c0   ::Vector{Float64}
    cr   ::Vector{Float64}
    phi  ::Vector{Float64}
    Hp   ::Vector{Float64}
    E    ::Vector{Float64}
    nu   ::Vector{Float64}
    rho0 ::Vector{Float64}
    n0   ::Vector{Float64}
    T    ::Vector{Float64}
    cp   ::Vector{Float64}
    k    ::Vector{Float64}
end

"""
    PointData(xp, ni; material=get_material()) -> PointData

Points at `xp` (D×nmp, `ni` per element and axis), at rest, all of one `material`.
"""
function PointData(xp::AbstractMatrix, ni::Integer; material=get_material())
    n = size(xp, 2)
    f(k) = fill(Float64(getfield(material, k)), n)
    return PointData(Float64.(xp), zeros(size(xp)), Int(ni), zeros(Int, n), f(:c0), f(:cr), f(:phi), f(:Hp), f(:E), f(:nu),
                     f(:rho0), f(:n0), f(:T), f(:cp), f(:k))
end
