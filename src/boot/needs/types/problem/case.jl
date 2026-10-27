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

"""
    inpolygon(x, z, v) -> Bool

Whether the point `(x, z)` lies inside the polygon with vertices `v` (`[x, z]` pairs, not repeated
at the end), by the even–odd rule, so concave polygons work.
"""
function inpolygon(x::Real, z::Real, v)
    inside, j = false, length(v)
    for i ∈ eachindex(v)
        xi, zi, xj, zj = v[i][1], v[i][2], v[j][1], v[j][2]
        if (zi > z) != (zj > z) && x < (xj - xi)*(z - zi)/(zj - zi) + xi
            inside = !inside
        end
        j = i
    end
    return inside
end

abstract type AbstractShape end
"""
    Polygon(vertices)
    Polygon(rings)
    Polygon(vertices, sizes)

A polygon in the x–z plane (`[x, z]` vertices, not repeated at the end; even–odd rule); in 3D it is
extruded along y. Several outlines (`rings`, each a vector of vertices) make one zone of several parts or
with holes, by the even–odd rule over all of them (inside an odd number of outlines: inside), as the drawing
window leaves a zone cut by the ones drawn after it. `vertices` holds every outline's vertices, one outline
after the other, and `sizes` how many each has (`polygon_rings`).
"""
struct Polygon <: AbstractShape
    vertices::Vector{Vector{Float64}}
    sizes::Vector{Int}
    function Polygon(v, sizes=nothing)
        nested = isnothing(sizes) && !isempty(v) && first(v) isa AbstractVector && !isempty(first(v)) && first(first(v)) isa AbstractVector
        rings  = nested ? collect(v) : [v]
        vs = [Float64.(collect(p)) for r ∈ rings for p ∈ r]
        ns = isnothing(sizes) ? [length(r) for r ∈ rings] : collect(Int, sizes)
        sum(ns) == length(vs) && all(≥(3), ns) && all(length.(vs) .== 2) ||
            throw(ArgumentError("a Polygon needs at least 3 [x, z] vertices per outline"))
        return new(vs, ns)
    end
end
"The outlines of a polygon, each a vector of its vertices."
function polygon_rings(s::Polygon)
    ends = cumsum(s.sizes)
    return [s.vertices[e-n+1:e] for (n, e) ∈ zip(s.sizes, ends)]
end
"""
    Block(lo, hi)

The axis-aligned box `lo ≤ x ≤ hi` (bounds included): `D` coordinates, or two (x, z) for a box
extruded along y in 3D.
"""
struct Block <: AbstractShape
    lo::Vector{Float64}
    hi::Vector{Float64}
    function Block(lo, hi)
        length(lo) == length(hi) || throw(ArgumentError("Block bounds need the same length"))
        all(lo .≤ hi) || throw(ArgumentError("Block needs lo ≤ hi"))
        return new(Float64.(collect(lo)), Float64.(collect(hi)))
    end
end
"""
    Ball(centre, radius)

The points at most `radius` from `centre` (bounds included): a disk in 2D, a sphere in 3D, or, with
an (x, z) centre in 3D, a cylinder along y.
"""
struct Ball <: AbstractShape
    centre::Vector{Float64}
    radius::Float64
    Ball(c, r) = r ≥ 0 ? new(Float64.(collect(c)), Float64(r)) : throw(ArgumentError("Ball needs radius ≥ 0"))
end
"""
    Region(f)

The points `x` (a vector of `D` coordinates) for which `f(x)` is true: any shape, from code. A
`Region` cannot be written to a TOML file.
"""
struct Region{F} <: AbstractShape
    f::F
end

# whether the point x (D coordinates) is in the shape; two-coordinate shapes act on (x, z)
xz(x) = (x[1], x[end])
function inside(s::Polygon, x)                                # even–odd over all its outlines
    odd, o = false, 0
    for n ∈ s.sizes
        inpolygon(x[1], x[end], view(s.vertices, o+1:o+n)) && (odd = !odd)
        o += n
    end
    return odd
end
function inside(s::Block, x)
    c = length(s.lo) == length(x) ? x : xz(x)
    length(c) == length(s.lo) || throw(ArgumentError("a Block of $(length(s.lo)) coordinates in $(length(x))D"))
    return all(i -> s.lo[i] ≤ c[i] ≤ s.hi[i], eachindex(s.lo))
end
function inside(s::Ball, x)
    c = length(s.centre) == length(x) ? x : xz(x)
    length(c) == length(s.centre) || throw(ArgumentError("a Ball of $(length(s.centre)) coordinates in $(length(x))D"))
    d2 = 0.0
    for i ∈ eachindex(s.centre)
        Δ = c[i] - s.centre[i]; d2 += Δ*Δ
    end
    return sqrt(d2) ≤ s.radius
end
inside(s::Region, x) = s.f(x)::Bool

"""
    Zone(name, shape; props...)
    Zone(name, shape, props::AbstractDict)

A part of the material: the material points whose sub-cell centre lies in `shape`, with `props`
replacing the case's material (`MATERIAL` keys, SI, angles in radians) and an initial velocity `v`.
A point belongs to the last zone that contains it.
"""
struct Zone
    name::String
    shape::AbstractShape
    props::Dict{Symbol,Any}
    function Zone(name::AbstractString, shape::AbstractShape, props::AbstractDict)
        p = Dict{Symbol,Any}(Symbol(k) => v for (k, v) ∈ props)
        bad = setdiff(keys(p), ZONE_PROPS)
        isempty(bad) || throw(ArgumentError("zone \"$name\": unknown properties $(collect(bad)); known: $(ZONE_PROPS)"))
        for (k, v) ∈ p
            p[k] = k == :v ? Float64.(collect(v)) : Float64(v)
        end
        return new(String(name), shape, p)
    end
end
Zone(name::AbstractString, shape::AbstractShape; props...) = Zone(name, shape, Dict{Symbol,Any}(props))

"""
    Case(; L, nel, ni=2, zones, material=(;), bcs, solver=(;), thermal=false, run=(;), width=nothing, nely=nothing)

A simulation case, the description a problem is built from: the domain `L` [m] with `nel` elements (2D or 3D; the lower corner at the origin)
and `ni`×`ni` material points per element; the `zones` holding material (`Zone`); the default
`material` (overrides of `MATERIAL`); the boundary condition of each side, `bcs` (one `[lower upper]`
row per axis, `:free`, `:roller` or `:fixed`; default roller; extra rows are ignored, as the solver
does in 2D); `solver` overrides of `get_default()`
(merged per field, `merge_config`; `solution` follows the workflow unless set); `thermal` (build the
thermal phase); `run` settings (`RUN`: `folder`, `te`, `tg`, `tep`, `workflow`, `lstps`). A 2D (x–z)
case with `width` and `nely` is extruded along y into 3D, its sides along y roller.
`setup_simulation` builds it, `save_case`/`load_case` store it as TOML.
"""
struct Case
    L        :: Vector{Float64}
    nel      :: Vector{Int}
    ni       :: Int
    width    :: Union{Nothing,Float64}
    nely     :: Union{Nothing,Int}
    zones    :: Vector{Zone}
    material :: NamedTuple
    bcs      :: Matrix{Symbol}
    solver   :: NamedTuple
    thermal  :: Bool
    run      :: NamedTuple
end
function Case(; L, nel, ni::Integer=2, zones=Zone[], material=(;), bcs=nothing, solver=(;), thermal::Bool=false,
                run=(;), width=nothing, nely=nothing)
    L, nel = Float64.(collect(L)), Int.(collect(nel))
    D = length(L)
    length(nel) == D && D ∈ (2, 3) || throw(ArgumentError("L and nel need 2 or 3 entries each"))
    all(L .> 0) && all(nel .≥ 1) && ni ≥ 1 || throw(ArgumentError("L > 0, nel ≥ 1 and ni ≥ 1 are needed"))
    isnothing(width) == isnothing(nely) || throw(ArgumentError("extrude to 3D with both width and nely"))
    isnothing(width) || D == 2 || throw(ArgumentError("only a 2D (x–z) case can be extruded along y"))
    isnothing(width) || (width > 0 && nely ≥ 1) || throw(ArgumentError("width > 0 and nely ≥ 1 are needed"))
    nt(x) = x isa AbstractDict ? (; (Symbol(k) => nt(v) for (k, v) ∈ x)...) : x
    m   = nt(material)
    bad = setdiff(keys(m), keys(MATERIAL))
    isempty(bad) || throw(ArgumentError("unknown material constants $(collect(bad)); known: $(keys(MATERIAL))"))
    mat = (; (k => Float64(get(m, k, getfield(MATERIAL, k))) for k ∈ keys(MATERIAL))...)
    s   = nt(solver)
    s isa NamedTuple || throw(ArgumentError("solver overrides are a NamedTuple or Dict of sections"))
    if haskey(s, :bcs)                                       # boundary conditions passed with the solver options
        isnothing(bcs) && (bcs = s.bcs.dirichlet)
        s = Base.structdiff(s, NamedTuple{(:bcs,)})
    end
    b   = isnothing(bcs) ? fill(:roller, D, 2) : Symbol.(bcs)
    # more rows than axes (the solver default has three): the first D, or x and z of an x–y–z matrix when extruded
    size(b, 1) > D && size(b, 2) == 2 && (b = !isnothing(width) && size(b, 1) == 3 ? b[[1, 3], :] : b[1:D, :])
    size(b) == (D, 2) && all(∈(BOUNDARY_KINDS), b) ||
        throw(ArgumentError("bcs is one [lower upper] row per axis ($(D)×2), each one of $(BOUNDARY_KINDS)"))
    merge_config(get_default(), (; (k => v for (k, v) ∈ pairs(s) if haskey(get_default(), k))...))   # unknown fields: error now
    return Case(L, nel, Int(ni), isnothing(width) ? nothing : Float64(width), isnothing(nely) ? nothing : Int(nely),
                       collect(Zone, zones), mat, b, s, thermal, case_run(run))
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

# the domain the case is built on (extruded along y when it has a width)
domain(s::Case) = isnothing(s.width) ? (s.L, s.nel) : ([s.L[1], s.width, s.L[2]], [s.nel[1], s.nely, s.nel[2]])

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
