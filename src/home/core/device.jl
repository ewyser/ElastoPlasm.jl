"""
    to_device(x, solver) -> x on the solver's device

`x` (a `Point`, `Mesh`, `Basis`, or a tuple of them) with its arrays on the solver's backend (`exec_backend`): the
same object, no copy, on the host (`CPU()`), otherwise `adapt(backend, x)`. The explicit workflows step on what it
returns; the host problem stays the one a `Simulation` holds, saves and plots (`to_host!`).
"""
function to_device(x, solver::AbstractSolver)
    backend = exec_backend(solver)
    return backend isa CPU ? x : Adapt.adapt(backend, x)
end

"""
    to_host!(host, dev) -> host

Copy `dev`'s arrays (`to_device`'s copy) back into `host`'s, field by field through the containers; nothing to do
when they are the same objects (the host). Plain-bits fields (scalars, static vectors, `MeshProperties`, `nothing`)
are equal on both sides, so they are skipped; an `SVector` is checked as one before as an array, which it also is.
"""
function to_host!(host, dev)
    host === dev && return host
    for f ∈ fieldnames(typeof(host))
        h, d = getfield(host, f), getfield(dev, f)
        if isbits(h)
            continue
        elseif h isa AbstractArray
            copyto!(h, d)
        else
            to_host!(h, d)
        end
    end
    return host
end

"""
    host_only(solver, name)

Refuse to run `name` (a workflow not yet ported to devices: dynamic relaxation, thermal) with a solver on a device,
with an `ArgumentError` naming the device, rather than failing deep inside a kernel.
"""
function host_only(solver::AbstractSolver, name::AbstractString)
    exec_backend(solver) isa CPU && return nothing
    throw(ArgumentError("$name runs on the host only for now, and this solver runs on " *
                        "$(first(values(solver.backend.exec)).name): set backend = (; select = \"host\")"))
end
