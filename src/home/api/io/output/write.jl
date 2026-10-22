export record!

"""
    Recorder

An open HDF5 output: the file, the recorded fields (`(field, object, dataset)`) and the element type
they are stored in. Made by `open_output`, written by `record!`, closed by `close`.
"""
struct Recorder
    file   :: HDF5.File
    fields :: Vector{Tuple{Function,String,HDF5.Dataset}}
    T      :: DataType
end
Base.close(r::Recorder) = isopen(r.file) && close(r.file)
nframes(r::Recorder) = length(r.file["frames/t"])

# a dataset of frames of shape `dims`, growing along its last axis, one chunk per frame
frames_dataset(f, path, T, dims) = create_dataset(f, path, T, ((dims..., 0), (dims..., -1)); chunk=(dims..., 1), deflate=3)
# set frame k of a dataset (the last axis), growing it to k frames
function set_frame!(d::HDF5.Dataset, k, values)
    HDF5.set_extent_dims(d, (size(d)[1:end-1]..., k))
    if ndims(d) == 1
        d[k] = values
    elseif ndims(d) == 2
        d[:, k] = values
    else
        d[:, :, k] = values
    end
end

"""
    open_output(file, problem; append=false, case_text=nothing) -> Recorder

The HDF5 output of a run: a new file, or with `append` an existing one, continued.
"""
open_output(file, p::Problem; append::Bool=false, case_text=nothing) =
    append && isfile(file) ? reopen_output(file, p) : create_output(file, p; case_text)

# a new output: the grid, the points as built, what was run, a dataset per field, and the first frame
function create_output(file, p::Problem; case_text=nothing)
    T = p.solver.output.single ? Float32 : Float64
    m, g = p.mpts, p.mesh
    mkpath(dirname(abspath(file)))
    f = h5open(file, "w")
    f["mesh/x"]   = T.(columns(g.x))
    f["mesh/L"]   = Float64.(collect(g.prprt.L))
    f["mesh/nel"] = Int.(collect(g.prprt.nel))
    f["mesh/h"]   = Float64.(collect(g.prprt.h))
    f["point/initial/Omega"] = T.(m.Ω₀)
    f["point/initial/n0"]    = T.(m.n₀)
    f["point/initial/rho0"]  = T.(m.s.ρ₀)
    a = attributes(f)
    a["ElastoPlasm"] = string(get_version())
    a["D"], a["nmp"], a["nno"] = length(g.prprt.L), Int(m.nmp), Int(g.prprt.nno[end])
    a["solution"], a["basis"], a["transfer"] = p.solver.solution, p.solver.basis.which, p.solver.basis.trsfr
    a["elastic"], a["plastic"] = p.solver.material.elastic, p.solver.material.plastic
    isnothing(case_text) || (a["case"] = case_text)
    frames_dataset(f, "frames/t", Float64, ())
    frames_dataset(f, "frames/step", Int, ())
    frames_dataset(f, "frames/loadstep", UInt8, ())                    # 1 when t is a load step (0/1: readable anywhere)
    frames_dataset(f, "frames/workflow", datatype(String), ())
    fields = [(field, object, frames_dataset(f, path, T, size(field(object == "point" ? m : g))))
              for (path, field, object) ∈ output_fields(p)]
    r = Recorder(f, fields, T)
    write_frame!(r, p.mpts, p.mesh, p.time.t[1], 0, false, "initial")
    return r
end

# an existing output, to continue: the fields it already records
function reopen_output(file, p::Problem)
    f = h5open(file, "r+")
    selected = output_fields(p)
    missing_ = [path for (path, _, _) ∈ selected if !haskey(f, path)]
    isempty(missing_) || @warn "Output: $(file) does not record $(missing_); not added when continuing it"
    fields = [(field, object, f[path]) for (path, field, object) ∈ selected if haskey(f, path)]
    T = isempty(fields) ? (p.solver.output.single ? Float32 : Float64) : eltype(fields[1][3])
    return Recorder(f, fields, T)
end

function write_frame!(r::Recorder, mpts, mesh, t, step, loadstep::Bool, workflow::AbstractString)
    k = nframes(r) + 1
    for (field, object, d) ∈ r.fields
        set_frame!(d, k, r.T.(field(object == "point" ? mpts : mesh)))
    end
    set_frame!(r.file["frames/t"], k, Float64(t))
    set_frame!(r.file["frames/step"], k, Int(step))
    set_frame!(r.file["frames/loadstep"], k, UInt8(loadstep))
    set_frame!(r.file["frames/workflow"], k, workflow)
end

"""
    record!(mpts, mesh, time, solver; step=0, loadstep=false)

Write a frame of the run's HDF5 output, if `execute!` has one open (with `solver.output.status`):
the selected fields at `time.t[1]` (a time, or with `loadstep` a load step) and the step count
`step`. Workflows call it at their checkpoints; it does nothing without an open output, and does not
write the same time twice in a row.
"""
function record!(mpts, mesh, time, solver; step::Integer=0, loadstep::Bool=false)
    solver.output.status || return nothing
    r = get(task_local_storage(), :output_recorder, nothing)
    isnothing(r) && return nothing
    t = Float64(time.t[1])
    n = nframes(r)
    n > 0 && r.file["frames/t"][n] == t && return nothing
    write_frame!(r, mpts, mesh, t, step, loadstep, get(task_local_storage(), :output_workflow, "?"))
    return nothing
end
