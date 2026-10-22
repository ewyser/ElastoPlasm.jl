# whether the problem has the phase `phase` of `object` ("point" or "mesh")
function has_phase(p::Problem, object, phase)
    phase ∈ ("base", "solid") && return true
    phase == "thermal"        && return !isnothing(object == "point" ? p.mpts.t : p.mesh.t)
    phase == "fluid"          && return object == "point" && !isnothing(p.mpts.f)
    return false
end

"""
    output_fields(problem) -> Vector{(path, entry, object)}

The fields `problem.solver.output.what` selects, as `("point/solid/p", entry, "point")` with the entry of
`quantities()`: unknown names are an error, fields of a phase the problem does not have are skipped with a
warning.
"""
function output_fields(p::Problem)
    selected, skipped = [], String[]
    for (object, phases) ∈ pairs(p.solver.output.what), (phase, names) ∈ pairs(phases), name ∈ names
        object, phase = string(object), string(phase)
        path = "$object/$phase/$name"
        e = quantity(path)
        isnothing(e) && throw(ArgumentError("unknown output field $path; $object/$phase has: " *
                                            "$(collect(String.(keys(something(quantity("$object/$phase"), (;))))))"))
        has_phase(p, object, phase) ? push!(selected, (path, e, object)) : push!(skipped, path)
    end
    isempty(skipped) || @warn "Output: skipped fields of phases this problem does not have: $skipped"
    return selected
end
