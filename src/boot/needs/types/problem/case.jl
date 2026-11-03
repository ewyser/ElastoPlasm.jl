# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# What ElastoPlasm reads from a case (EasyDraw's `Case`: a layout with its settings): its tables (the material, the zone
# properties, the run settings, the boundary kinds, the workflows), its physics (`case_physics`), its settings written
# in ElastoPlasm's terms (`case_settings`), and the material points a problem starts from (`PointData`)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractShape, Polygon, Block, Ball, Region, Zone, Layout, Settings, load_layout, load_settings, save_settings,
       Case, save_case, load_case, case_settings, case_physics, PointData, get_material, inpolygon

"""
    get_material() -> NamedTuple

The default material of a case, in SI units with angles in radians: Young's modulus `E` [Pa], Poisson's ratio `nu`,
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
# the shapes (`Polygon`, `Block`, `Ball`, `Region`, their `inside` rule), `Zone`, `Layout`, `Settings` and `Case` are
# EasyDraw's (lib/EasyDraw), imported in boot.jl; a case's zones may set only `get_zone_properties()`, checked by
# `case_physics`

"""
    case_settings(; properties=Dict(), material=(;), bcs=nothing, solver=(;), thermal=false, run=(;)) -> Settings

A case's settings written in ElastoPlasm's terms, for `Case(layout, case_settings(; …))` (EasyDraw's `Case`: the
layout is the geometry, these are its values):
- `properties`: each zone's constants, by zone name: `Dict("rock" => (; E=5e7, phi=30π/180))`, any `get_material()`
  key or an initial velocity `v`;
- `material`: overrides of `get_material()` (SI, angles in radians), the settings' `defaults`;
- `bcs`: the boundary conditions, one `[lower upper]` row per axis (x, z; or x, y, z), each `:free`, `:roller` or
  `:fixed`, or `solver.bcs.dirichlet` as the solver options give them; the settings' `sides`, a roller side left out
  (a side not set is roller);
- `solver`: overrides of `get_default()`, the physics (`material`, `grf`) and the numerics alike, written in their
  file form (`option_labels`);
- `thermal`, `run`: the run settings, every `get_run()` one (`case_run`), with `thermal` when true.

# Throws
- `ArgumentError`: for an unknown or invalid option or run setting, or a `bcs` that is not 2×2 or 3×2.
"""
function case_settings(; properties=Dict(), material=(;), bcs=nothing, solver=(;), thermal::Bool=false, run=(;))
    s = named_tuple(solver)
    if haskey(s, :bcs)                                       # boundary conditions passed with the solver options
        isnothing(bcs) && (bcs = s.bcs.dirichlet)
        s = Base.structdiff(s, NamedTuple{(:bcs,)})
    end
    sides    = isnothing(bcs) ? Dict{String,String}() : case_sides(Symbol.(bcs))
    defaults = Dict{String,Any}(string(k) => v for (k, v) ∈ pairs(named_tuple(material)))
    options  = option_labels(config_dict(s))
    case_options(options)                                    # checked: an unknown option errors here, not at setup
    r = config_dict(case_run(run))
    thermal && (r["thermal"] = true)
    return EasyDraw.Settings(; zones=properties, defaults, sides, options, run=r)
end
# the sides of a `bcs` matrix, one [lower upper] row per axis (its rows give the dimension), a roller side left out:
# what a side not set is (`case_bcs`)
function case_sides(bcs::AbstractMatrix)
    size(bcs, 1) ∈ (2, 3) && size(bcs, 2) == 2 || throw(ArgumentError("bcs is one [lower upper] row per axis (2×2 or 3×2)"))
    return Dict(side => string(bcs[axis, k]) for (side, (axis, k)) ∈ side_axes(size(bcs, 1)) if bcs[axis, k] != :roller)
end

"""
    case_physics(case) -> NamedTuple

What ElastoPlasm reads from a case's settings as its physics, checked: what `Problem(case, solver)` builds the points
and the mesh from. A case (EasyDraw's `Case`: a layout with its settings) gives its values no meaning; this does.

# Fields
- `material::NamedTuple`: the material constants: `get_material()`, with the settings' `defaults`.
- `properties::Dict{String,NamedTuple}`: each zone's own constants, by zone name: `get_material()` keys, or an
  initial velocity `v`; where a zone does not set one, the material's.
- `bcs::Matrix{Symbol}`: the boundary condition of each side, one `[lower upper]` row per axis (x, then y in 3D,
  then z): the settings' `sides`, `"left"`/`"right"`, `"front"`/`"back"`, `"bottom"`/`"top"`, each `:free`,
  `:roller` (a side not set) or `:fixed`.
- `laws::NamedTuple`: the constitutive laws, `(; elastic, plastic)`: the options `material.elastic` (`"hencky"`,
  `"improved_hencky"` or `"hypoelastic"`) and `material.plastic` (`"DP"` or `"VM"`), `get_default()`'s where not set.
- `grf::NamedTuple`: the random field of the cohesion, `(; status, covariance, param)`: the options `grf`,
  `get_default()`'s where not set.
- `thermal::Bool`: whether to build the thermal phase: the run setting `thermal`, false when absent.

# Throws
- `ArgumentError`: when a zone sets something other than a material constant or a velocity `v`, a default is not a
  material constant, a side's kind is not one of `get_boundary_kinds()`, or an option or a run setting is unknown or
  invalid (all are checked, the numerical ones too, which `Solver(case)` reads).
"""
function case_physics(case::Case)
    l, s = case.layout, case.settings
    options = case_options(s.options)                        # every option, checked
    case_run(l, s)                                           # every run setting, checked
    d = get_default()
    laws = haskey(options, :material) ? merge_config(d.material, options.material; path="material.") : d.material
    grf  = haskey(options, :grf) ? merge_config(d.grf, options.grf; path="grf.") : d.grf
    return (; material=case_material(s.defaults), properties=case_properties(s.zones), bcs=case_bcs(s.sides, length(l.L)),
              laws, grf, thermal=Bool(get(s.run, "thermal", false)))
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
