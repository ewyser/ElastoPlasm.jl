export plot_output, animate_output, plot_series

# how a recorded field is shown: the live plots' settings (`get_mpts_variable_config`) where one matches
const OUTPUT_PLOT_NAMES = Dict("p" => "P", "epII" => "epII", "J" => "J", "n" => "n", "umag" => "u", "vmag" => "v", "T" => "T")
function output_style(field::AbstractString, component)
    name = last(split(field, "/"))
    cfg  = get_mpts_variable_config()
    key  = get(OUTPUT_PLOT_NAMES, name, nothing)
    isnothing(component) && !isnothing(key) && haskey(cfg, key) && return cfg[key]
    sub  = isnothing(component) ? "" : "[$component]"
    return (; label=name * sub, unit="", scale=1.0, name=name * sub, cb=:viridis, cblim=nothing)
end
# one frame's values of a field: a scalar per point, a component, or the magnitude of a vector
function frame_values(v::AbstractArray, component)
    ndims(v) == 1 && return v
    isnothing(component) || return v[component, :]
    size(v, 1) ≤ 3 && return vec(sqrt.(sum(abs2, v; dims=1)))
    throw(ArgumentError("a field of $(size(v, 1)) components: pass component=…"))
end
frame_title(name, t, loadstep) = loadstep ? "$(uppercasefirst(name)), at load step $(round(Int, t))" :
                                            "$(uppercasefirst(name)), at " * L" t = " * string(round(t; digits=1)) * " [s]"

"""
    plot_output(file, field="point/solid/p"; frame=lastindex, component=nothing, dpi=500, file_out=nothing) -> Plot

Draw one recorded frame of a point field from an HDF5 output (`execute!` with `output.status`),
as the live plots do: the points at their recorded positions (`point/base/x`; in 3D the x–z plane),
coloured by the field, with the same colour map, unit and limits where the field has a live-plot
setting, titled "at t = …" or "at load step …". A vector field shows its magnitude or `component`.
Saves the figure when `file_out` is given. Needs no `Problem`: only the file.
"""
function plot_output(file::AbstractString, field::AbstractString="point/solid/p"; frame=nothing, component=nothing,
                     dpi::Integer=500, file_out=nothing)
    startswith(field, "point/") || throw(ArgumentError("plot_output draws point fields; $(field) is not one"))
    o  = load_frame_plotdata(file, field, frame)
    st = output_style(field, component)
    d  = frame_values(o.v, component) .* st.scale
    clim = isnothing(st.cblim) ? (minimum(d) == maximum(d) ? (minimum(d) - 1, maximum(d) + 1) : extrema(d)) : st.cblim
    xz   = (o.x[1, :], o.x[end, :])
    Lx, Lz = o.L[1], o.L[end]
    dims = (dpi, round(Int, dpi*Lz/Lx))
    ms   = min(dims[1]/Lx, dims[2]/Lz) * o.spacing
    config_plot()
    p = plot(xz...; seriestype=:scatter, marker_z=d, markersize=ms, markershape=:circle, markerstrokewidth=0.75,
             xlabel=L"$x-$direction"*" [m]", ylabel=L"$z-$direction"*" [m]", label=false, legend=false, colorbar=true,
             colorbar_title="$(st.label) $(st.unit)", color=st.cb, clim=clim, xlim=(0, Lx), ylim=(0, Lz),
             title=frame_title(st.name, o.t, o.loadstep), aspect_ratio=1, size=dims)
    isnothing(file_out) || savefig(p, file_out)
    return p
end
# the positions, values and frame data `plot_output` needs, read from the file
function load_frame_plotdata(file, field, frame)
    h5open(file, "r") do f
        haskey(f, "point/base/x") || throw(ArgumentError("$(file) has no point positions (record point.base x)"))
        haskey(f, field) || throw(ArgumentError("$(file) has no field $(field)"))
        n = length(f["frames/t"]); k = something(frame, n)
        1 ≤ k ≤ n || throw(ArgumentError("frame $k of $n"))
        d = f[field]
        D = Int(read(attributes(f)["D"]))
        return (; x=f["point/base/x"][:, :, k], v=(ndims(d) == 2 ? d[:, k] : d[:, :, k]), t=f["frames/t"][k],
                  loadstep=f["frames/loadstep"][k] == 1, L=read(f["mesh/L"]), spacing=(sum(read(f["point/initial/Omega"]))/length(f["point/initial/Omega"]))^(1/D))
    end
end

"""
    animate_output(file, field="point/solid/p"; fps=10, component=nothing, dpi=300, file_out) -> String

A GIF (or MP4, from `file_out`'s extension) of every recorded frame of a point field, each drawn as
`plot_output` does; returns the written path (default `<output>_<field>.gif` beside the file).
"""
function animate_output(file::AbstractString, field::AbstractString="point/solid/p"; fps::Integer=10, component=nothing,
                        dpi::Integer=300, file_out=replace(file, r"\.h5$" => "") * "_" * replace(field, "/" => "_") * ".gif")
    n = h5open(f -> length(f["frames/t"]), file, "r")
    anim = @animate for k ∈ 1:n
        plot_output(file, field; frame=k, component, dpi)
    end
    endswith(lowercase(file_out), ".mp4") ? mp4(anim, file_out; fps) : gif(anim, file_out; fps)
    return file_out
end

"""
    plot_series(file, field="point/solid/p"; point, component=nothing, file_out=nothing) -> Plot

The history of a recorded field at one point (a node for a mesh field) against time, or the load step.
"""
function plot_series(file::AbstractString, field::AbstractString="point/solid/p"; point::Integer, component=nothing, file_out=nothing)
    s  = output_series(file, field; point, component)
    st = output_style(field, component)
    v  = s.values isa AbstractMatrix ? (size(s.values, 1) ≤ 3 ? vec(sqrt.(sum(abs2, s.values; dims=1))) : throw(ArgumentError("pass component=…"))) : s.values
    config_plot()
    p = plot(s.t, v .* st.scale; marker=:circle, label=false, xlabel=s.loadstep ? "load step" : "t [s]",
             ylabel="$(st.label) $(st.unit)", title="$(uppercasefirst(st.name)) at $(field[1:5] == "point" ? "point" : "node") $(point)", size=(500, 300))
    isnothing(file_out) || savefig(p, file_out)
    return p
end
