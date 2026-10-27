export output_series, output_frame, load_output

# the paths of the recorded fields, "point/solid/p" etc.
function recorded_fields(f)
    paths = String[]
    for object ∈ ("point", "mesh")
        haskey(f, object) || continue
        for phase ∈ keys(f[object])
            (phase == "initial" || !(f[object][phase] isa HDF5.Group)) && continue
            append!(paths, ["$object/$phase/$name" for name ∈ keys(f[object][phase])])
        end
    end
    return paths
end
# frame k of a dataset (the last axis)
frame(d, k) = ndims(d) == 1 ? d[k] : ndims(d) == 2 ? d[:, k] : d[:, :, k]
# what a recorded field is, from its attributes: readable name, SI unit, LaTeX label (an older file: its key)
function field_info(f, path)
    a = attributes(f[path]); get_attr(k, d) = haskey(a, k) ? read(a[k]) : d
    return (; name=get_attr("name", basename(path)), unit=get_attr("unit", ""), label=get_attr("label", nothing))
end

"""
    output_series(file, field; point=nothing, component=nothing) -> (; t, values, loadstep)

The history of a recorded field (`"point/solid/p"`, `"point/base/x"`, `"mesh/solid/v"`, …) against
`t` (a time, or a load step when `loadstep`): at one `point` (a node for a mesh field), or every point
without one; one `component` of a vector or deviator field, or all of them without one.
"""
function output_series(file::AbstractString, field::AbstractString; point=nothing, component=nothing)
    h5open(file, "r") do f
        haskey(f, field) || throw(ArgumentError("$(file) has no field $(field); recorded: $(recorded_fields(f))"))
        d = f[field]
        if ndims(d) == 2                                        # a scalar: points × frames
            isnothing(component) || throw(ArgumentError("$(field) has no components"))
            values = isnothing(point) ? read(d) : d[point, :]
        else                                                    # components × points × frames
            c = something(component, Colon())
            values = isnothing(point) ? d[c, :, :] : d[c, point, :]
        end
        loadstep = read(f["frames/loadstep"]) .== 1
        return (; t=read(f["frames/t"]), values, loadstep=length(loadstep) > 1 && all(loadstep[2:end]))
    end
end

"""
    output_frame(file, k) -> (; t, step, loadstep, workflow, fields, info)

Every recorded field at frame `k` (`fields`: path → that frame; `info`: path → its name, SI unit, label).
"""
function output_frame(file::AbstractString, k::Integer)
    h5open(file, "r") do f
        n = length(f["frames/t"])
        1 ≤ k ≤ n || throw(ArgumentError("frame $k of $n"))
        fields = Dict{String,Any}(path => frame(f[path], k) for path ∈ recorded_fields(f))
        info   = Dict{String,Any}(path => field_info(f, path) for path ∈ recorded_fields(f))
        return (; t=f["frames/t"][k], step=f["frames/step"][k], loadstep=f["frames/loadstep"][k] == 1,
                  workflow=f["frames/workflow"][k], fields, info)
    end
end

"""
    load_output(file) -> (; frames, fields, info, mesh, initial, attributes)

The whole output: `frames` (`t`, `step`, `loadstep`, `workflow` per frame), `fields` (path → array,
frames last), `info` (path → the field's name, SI unit and label), the grid `mesh`, the `initial` point
properties and the file's `attributes`.
"""
function load_output(file::AbstractString)
    h5open(file, "r") do f
        group(g) = Dict{String,Any}(k => read(g[k]) for k ∈ keys(g))
        frames  = (; t=read(f["frames/t"]), step=read(f["frames/step"]), loadstep=read(f["frames/loadstep"]) .== 1,
                     workflow=read(f["frames/workflow"]))
        fields  = Dict{String,Any}(path => read(f[path]) for path ∈ recorded_fields(f))
        info    = Dict{String,Any}(path => field_info(f, path) for path ∈ recorded_fields(f))
        mesh    = (; x=read(f["mesh/x"]), L=read(f["mesh/L"]), nel=read(f["mesh/nel"]), h=read(f["mesh/h"]))
        attrs   = Dict{String,Any}(k => read(attributes(f)[k]) for k ∈ keys(attributes(f)))
        return (; frames, fields, info, mesh, initial=group(f["point/initial"]), attributes=attrs)
    end
end
