export OUTPUT_FIELDS

columns(v) = reduce(hcat, v)                                 # Vector{SVector{D}} → D × n

"""
    OUTPUT_FIELDS

What an HDF5 output can record, shaped like `output.what`: object (`"point"`, `"mesh"`) → phase
(`"base"` for the point fields outside any phase, `"solid"`, `"fluid"`, `"thermal"`) → name → a
function giving one frame of the field (a vector over points or nodes, or components × count), raw SI.
To record a quantity of your own, add an entry: `OUTPUT_FIELDS["point"]["solid"]["J2"] = m -> …`.
"""
const OUTPUT_FIELDS = Dict{String,Dict{String,Dict{String,Function}}}(
    "point" => Dict(
        "base" => Dict(
            "x"     => m -> columns(m.x),
            "J"     => m -> m.J,
            "n"     => m -> m.n,
            "Omega" => m -> m.Ω,
            "ell"   => m -> columns(m.ℓ),
        ),
        "solid" => Dict(
            "u"     => m -> columns(m.s.u),
            "v"     => m -> columns(m.s.v),
            "p"     => m -> getfield.(m.s.σᵢ, :p),
            "dev"   => m -> columns([σ.dev for σ ∈ m.s.σᵢ]),
            "epII"  => m -> getindex.(m.s.ϵpII, 1),
            "epV"   => m -> m.s.ϵpV,
            "rho"   => m -> m.s.ρ,
            "q"     => m -> [sqrt(3*get_J2(σ)) for σ ∈ m.s.σᵢ],     # von Mises stress
            "umag"  => m -> [sqrt(sum(abs2, u)) for u ∈ m.s.u],
            "vmag"  => m -> [sqrt(sum(abs2, v)) for v ∈ m.s.v],
        ),
        "fluid"   => Dict{String,Function}(),                        # until the fluid phase has fields
        "thermal" => Dict("T" => m -> m.t.T),
    ),
    "mesh" => Dict(
        "solid" => Dict(
            "v" => g -> columns(g.s.v),
            "a" => g -> columns(g.s.a),
            "m" => g -> g.s.m,
        ),
        "thermal" => Dict("T" => g -> g.t.T),
    ),
)

# whether the problem has the phase `phase` of `object` ("point" or "mesh")
function has_phase(p::Problem, object, phase)
    phase ∈ ("base", "solid") && return true
    phase == "thermal"        && return !isnothing(object == "point" ? p.mpts.t : p.mesh.t)
    phase == "fluid"          && return object == "point" && !isnothing(p.mpts.f)
    return false
end

"""
    output_fields(problem) -> Vector{(path, field, object)}

The fields `problem.solver.output.what` selects, as `("point/solid/p", function, "point")`: unknown
names are an error, fields of a phase the problem does not have are skipped with a warning.
"""
function output_fields(p::Problem)
    selected, skipped = [], String[]
    for (object, phases) ∈ pairs(p.solver.output.what), (phase, names) ∈ pairs(phases), name ∈ names
        object, phase = string(object), string(phase)
        path = "$object/$phase/$name"
        known = get(get(OUTPUT_FIELDS, object, Dict()), phase, Dict())
        haskey(known, name) || throw(ArgumentError("unknown output field $path; $object/$phase has: $(sort(collect(keys(known))))"))
        has_phase(p, object, phase) ? push!(selected, (path, known[name], object)) : push!(skipped, path)
    end
    isempty(skipped) || @warn "Output: skipped fields of phases this problem does not have: $skipped"
    return selected
end
