export OUTPUT_FIELDS, record!, output_series, output_frame, load_output

# ---- what can be recorded: (object, phase, name) → (the values of one frame, components)

vecs(v) = reduce(hcat, v)                                   # Vector{SVector{D}} → D × n
"""
    OUTPUT_FIELDS

The fields an HDF5 output may record (`output.what`), by object (`:point`, `:mesh`), phase (`:base`
for the point fields outside any phase, `:solid`, `:fluid`, `:thermal`) and name: a function of the
`mpts` (or `mesh`) giving one frame (a vector, or a matrix of components × points/nodes), raw SI.
Add an entry to record a quantity of your own.
"""
const OUTPUT_FIELDS = Dict{Tuple{Symbol,Symbol,String},Function}(
    (:point, :base,  "x")     => m -> vecs(m.x),
    (:point, :base,  "J")     => m -> m.J,
    (:point, :base,  "n")     => m -> m.n,
    (:point, :base,  "Omega") => m -> m.Ω,
    (:point, :base,  "ell")   => m -> vecs(m.ℓ),
    (:point, :solid, "u")     => m -> vecs(m.s.u),
    (:point, :solid, "v")     => m -> vecs(m.s.v),
    (:point, :solid, "p")     => m -> getfield.(m.s.σᵢ, :p),
    (:point, :solid, "dev")   => m -> vecs([σ.dev for σ ∈ m.s.σᵢ]),
    (:point, :solid, "epII")  => m -> getindex.(m.s.ϵpII, 1),
    (:point, :solid, "epV")   => m -> m.s.ϵpV,
    (:point, :solid, "rho")   => m -> m.s.ρ,
    (:point, :solid, "q")     => m -> [sqrt(3*get_J2(σ)) for σ ∈ m.s.σᵢ],      # von Mises stress
    (:point, :solid, "umag")  => m -> [sqrt(sum(abs2, u)) for u ∈ m.s.u],
    (:point, :solid, "vmag")  => m -> [sqrt(sum(abs2, v)) for v ∈ m.s.v],
    (:point, :thermal, "T")   => m -> m.t.T,
    (:mesh,  :solid, "v")     => g -> vecs(g.s.v),
    (:mesh,  :solid, "a")     => g -> vecs(g.s.a),
    (:mesh,  :solid, "m")     => g -> g.s.m,
    (:mesh,  :thermal, "T")   => g -> g.t.T,
)
# whether a problem has the phase a field belongs to
has_phase(mpts, mesh, obj, phase) = phase == :base || (obj == :point ? (phase == :solid || (phase == :thermal && !isnothing(mpts.t)) ||
                                                                         (phase == :fluid && !isnothing(mpts.f))) :
                                                                       (phase == :solid || (phase == :thermal && !isnothing(mesh.t))))

# ---- the recorder: one open file, one dataset per field, frames appended

"""
    Recorder

An open HDF5 output: the file, the recorded fields (`(path, extractor, object, dataset)`), the frame
datasets, the number of frames written and the element type. Made by `open_output`, written by
`record!`, closed by `close`.
"""
mutable struct Recorder
    file   :: HDF5.File
    fields :: Vector{Tuple{String,Function,Symbol,HDF5.Dataset}}
    frames :: NamedTuple
    n      :: Int
    T      :: DataType
    names  :: Vector{String}                                  # the workflows, by their frame code
end
Base.close(r::Recorder) = isopen(r.file) && close(r.file)

# a dataset of frames of shape `dims` growing along its last (frame) axis, chunked one frame at a time
grow_dataset(parent, name, T, dims::Tuple) = create_dataset(parent, name, T, ((dims..., 0), (dims..., -1)); chunk=(dims..., 1), deflate=3)
# write frame `k` of a dataset
function append_frame!(d::HDF5.Dataset, k::Int, values)
    s = size(d)
    HDF5.set_extent_dims(d, (s[1:end-1]..., k))
    ndims(d) == 1 ? (d[k] = values) : ndims(d) == 2 ? (d[:, k] = values) : (d[:, :, k] = values)
end

"""
    open_output(file, problem; append=false, case_text=nothing) -> Recorder

Open the HDF5 output of a problem for its `solver.output` selection: a new file (the mesh grid, the
initial point properties, attributes describing the run, and a frame of the current state), or with
`append=true` an existing one, whose fields and frames are continued. Fields of a phase the problem
does not have are skipped with a warning.
"""
function open_output(file::AbstractString, p::Problem; append::Bool=false, case_text=nothing)
    out  = p.solver.output
    T    = out.single ? Float32 : Float64
    sel  = [(Symbol(obj), Symbol(phase), String(name)) for (obj, phases) ∈ pairs(out.what) for (phase, names) ∈ pairs(phases) for name ∈ names]
    keep = filter(s -> has_phase(p.mpts, p.mesh, s[1], s[2]), sel)
    length(keep) < length(sel) && @warn "Output: skipped fields of phases this problem does not have: $(setdiff(sel, keep))"
    bad = filter(s -> !haskey(OUTPUT_FIELDS, s), keep)
    isempty(bad) || throw(ArgumentError("unknown output fields $(bad); known: $(sort(collect(keys(OUTPUT_FIELDS))))"))
    path(s) = "$(s[1])/$(s[2])/$(s[3])"
    if append && isfile(file)
        f = h5open(file, "r+")
        frames = (; t=f["frames/t"], step=f["frames/step"], kind=f["frames/kind"], workflow=f["frames/workflow"])
        fields = [(path(s), OUTPUT_FIELDS[s], s[1], f[path(s)]) for s ∈ keep if haskey(f, path(s))]
        length(fields) < length(keep) && @warn "Output: $(file) does not record $(setdiff(path.(keep), first.(fields))); not added on append"
        names = haskey(attributes(f["frames/workflow"]), "names") ? read(attributes(f["frames/workflow"])["names"]) : String[]
        return Recorder(f, fields, frames, length(frames.t), isempty(fields) ? T : eltype(fields[1][4]), names)
    end
    mkpath(dirname(abspath(file)))
    f = h5open(file, "w")
    m, g = p.mpts, p.mesh
    # once: the grid, the points as built, and what was run
    f["mesh/x"] = T.(vecs(g.x)); f["mesh/L"] = collect(Float64.(g.prprt.L)); f["mesh/nel"] = collect(Int.(g.prprt.nel)); f["mesh/h"] = collect(Float64.(g.prprt.h))
    f["point/initial/Omega"] = T.(m.Ω₀); f["point/initial/n0"] = T.(m.n₀); f["point/initial/rho0"] = T.(m.s.ρ₀)
    a = attributes(f)
    a["ElastoPlasm"] = string(get_version()); a["D"] = length(g.prprt.L); a["nmp"] = Int(m.nmp); a["nno"] = Int(g.prprt.nno[end])
    a["solution"] = p.solver.solution; a["basis"] = p.solver.basis.which; a["transfer"] = p.solver.basis.trsfr
    a["elastic"] = p.solver.material.elastic; a["plastic"] = p.solver.material.plastic
    isnothing(case_text) || (a["case"] = case_text)
    frames = (; t=grow_dataset(f, "frames/t", Float64, ()), step=grow_dataset(f, "frames/step", Int, ()),
                kind=grow_dataset(f, "frames/kind", Int, ()), workflow=grow_dataset(f, "frames/workflow", Int, ()))
    attributes(frames.kind)["meaning"] = "0: t is a time [s]; 1: t is a load step (dynamic relaxation)"
    fields = map(keep) do s
        v = OUTPUT_FIELDS[s](s[1] == :point ? m : g)
        (path(s), OUTPUT_FIELDS[s], s[1], grow_dataset(f, path(s), T, size(v)))
    end
    r = Recorder(f, fields, frames, 0, T, String[])
    write_frame!(r, m, g, p.time.t[1], 0, 0, "initial")
    return r
end

function write_frame!(r::Recorder, mpts, mesh, t, step, kind, workflow::AbstractString)
    k = r.n + 1
    w = findfirst(==(workflow), r.names)
    isnothing(w) && (push!(r.names, workflow); w = length(r.names))
    for (_, f, obj, d) ∈ r.fields
        append_frame!(d, k, r.T.(f(obj == :point ? mpts : mesh)))
    end
    append_frame!(r.frames.t, k, Float64(t)); append_frame!(r.frames.step, k, Int(step))
    append_frame!(r.frames.kind, k, Int(kind)); append_frame!(r.frames.workflow, k, w)
    a = attributes(r.frames.workflow)
    haskey(a, "names") && delete_attribute(r.frames.workflow, "names")
    a["names"] = r.names
    r.n = k
    return r
end

"""
    record!(mpts, mesh, time, solver; step=0, loadstep=false)

Write a frame of the run's HDF5 output, when `elastoplasm` has one open and `solver.output.status`:
the selected fields at time `time.t[1]` (or, with `loadstep=true`, at a load step), with the step
count `step`. Does nothing otherwise, so workflows call it at every checkpoint. A frame at the same
time as the previous one is not written twice.
"""
function record!(mpts, mesh, time, solver; step::Integer=0, loadstep::Bool=false)
    solver.output.status || return nothing
    r = get(task_local_storage(), :elastoplasm_output, nothing)
    isnothing(r) && return nothing
    t = Float64(time.t[1])
    r.n > 0 && r.frames.t[r.n] == t && return nothing
    write_frame!(r, mpts, mesh, t, step, loadstep ? 1 : 0, get(task_local_storage(), :elastoplasm_workflow, "?"))
    return nothing
end

"""
    checkpoints(time, solver, stop) -> Vector{Tuple{T, Bool, Bool}}

The stops of an explicit workflow from the current time to `stop`: `(t, plot, record)`, on the
`plot.freq` grid (the live plot) and, when `output.status`, the `output.every` grid (a frame).
Without output the stops are the plot grid alone, as before the output existed.
"""
function checkpoints(time::Time{T1,T2}, solver, stop) where {T1,T2}
    plots = T2.(sort(collect(time.t[1]:solver.plot.freq:stop)))
    solver.output.status || return [(T, true, false) for T ∈ plots]
    recs  = T2.(sort(collect(time.t[1]:solver.output.every:stop)))
    return [(T, T ∈ plots, T ∈ recs) for T ∈ sort(unique([plots; recs]))]
end

# ---- reading

"""
    output_series(file, field; point=nothing, component=nothing) -> (; t, values, loadstep)

The history of a recorded field (`"point/solid/p"`, `"point/base/x"`, `"mesh/solid/v"`, …) against the
frame axis `t` (a time, or a load step when `loadstep`): at one `point` (a node index for a mesh field)
a vector over frames, or with `component` (of a vector or deviator field) its component; without
`point`, every point (count × frames, or components × count × frames).
"""
function output_series(file::AbstractString, field::AbstractString; point=nothing, component=nothing)
    h5open(file, "r") do f
        haskey(f, field) || throw(ArgumentError("$(file) has no field $(field); recorded: $(recorded_fields(f))"))
        d, t, kind = f[field], read(f["frames/t"]), read(f["frames/kind"])
        v = if isnothing(point)
            isnothing(component) ? read(d) : (ndims(d) == 3 ? d[component, :, :] : throw(ArgumentError("$(field) has no components")))
        else
            ndims(d) == 2 ? d[point, :] : isnothing(component) ? d[:, point, :] : d[component, point, :]
        end
        return (; t, values=v, loadstep=!isempty(kind) && all(==(1), kind[2:end]))
    end
end

"""
    output_frame(file, k) -> (; t, step, loadstep, workflow, fields)

Every recorded field at frame `k` (`fields`: path → array of that frame).
"""
function output_frame(file::AbstractString, k::Integer)
    h5open(file, "r") do f
        n = length(f["frames/t"])
        1 ≤ k ≤ n || throw(ArgumentError("frame $k of $n"))
        names = read(attributes(f["frames/workflow"])["names"])
        fields = Dict{String,Array}(p => (d = f[p]; ndims(d) == 2 ? d[:, k] : d[:, :, k]) for p ∈ recorded_fields(f))
        return (; t=f["frames/t"][k], step=f["frames/step"][k], loadstep=f["frames/kind"][k] == 1,
                  workflow=names[f["frames/workflow"][k]], fields)
    end
end

"""
    load_output(file) -> NamedTuple

The whole output: `frames` (`t`, `step`, `kind`, `workflow` names per frame), `fields` (path → array,
frames last), `mesh` (the grid), `initial` (the point properties as built) and the file's `attributes`.
"""
function load_output(file::AbstractString)
    h5open(file, "r") do f
        names = read(attributes(f["frames/workflow"])["names"])
        frames = (; t=read(f["frames/t"]), step=read(f["frames/step"]), kind=read(f["frames/kind"]),
                    workflow=names[read(f["frames/workflow"])])
        fields  = Dict{String,Array}(p => read(f[p]) for p ∈ recorded_fields(f))
        mesh    = (; x=read(f["mesh/x"]), L=read(f["mesh/L"]), nel=read(f["mesh/nel"]), h=read(f["mesh/h"]))
        initial = Dict{String,Array}(k => read(f["point/initial/$k"]) for k ∈ keys(f["point/initial"]))
        attrs   = Dict{String,Any}(k => read(attributes(f)[k]) for k ∈ keys(attributes(f)))
        return (; frames, fields, mesh, initial, attributes=attrs)
    end
end
# the paths of the recorded fields: point/<phase>/<name>, mesh/<phase>/<name>
recorded_fields(f) = [join([obj, phase, name], "/") for obj ∈ ("point", "mesh") if haskey(f, obj)
                      for phase ∈ keys(f[obj]) if f[obj][phase] isa HDF5.Group && phase != "initial"
                      for name ∈ keys(f[obj][phase])]

