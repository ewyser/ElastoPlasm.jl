export get_solver

"""
    exec_backend(x) -> KernelAbstractions backend

The backend a solver runs on, for the solver or its configuration (`instr`): its first selected device's
(`backend.exec`, set by `get_solver` through UnifiedBackend), `CPU()` for the host. The kernels are built on it
(`init_ignite`, `init_mapsto`, `init_update`, `init_implicit`).
"""
function exec_backend(x)
    return first(values(x.backend.exec)).Backend
end

"""
    solver_config(solver) -> NamedTuple

A solver's configuration, the sections `build_solver` takes: what `save_simulation` keeps of a solver, without the
selected devices (`backend.exec`) or the kernels (`cairn`), which `load_simulation` rebuilds where the file is loaded.
`stored_solver(solver)` is the solver holding only it (no kernels: `cairn = nothing`).
"""
function solver_config(s::AbstractSolver)
    return (; solution=s.solution, dtype=s.dtype, basis=s.basis, stab=s.stab, bcs=s.bcs, grf=s.grf, material=s.material,
              nonloc=s.nonloc, plot=s.plot, output=s.output, backend=(; select=s.backend.select, distributed=s.backend.distributed))
end
function stored_solver(s::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    c = solver_config(s)
    return Base.typename(S).wrapper{T1,T2,D}(c.solution, c.dtype, c.basis, c.stab, c.bcs, c.grf, c.material, c.nonloc,
                                             c.plot, c.output, c.backend, nothing)
end

"""
    merge_config(ref, user; path="") -> NamedTuple

`user`'s settings over `ref`'s, section by section: nested sections are merged field by field
(`(; basis=(; which="gimpm"))` changes only `basis.which`), anything else (values, vectors such as
`plot.what`) replaces the default. A field `ref` does not have is an error naming the valid ones
(`backend.exec`, set by `get_solver`, is accepted).
"""
function merge_config(ref::NamedTuple, user::NamedTuple; path::String="")
    for k ∈ keys(user)
        (haskey(ref, k) || (path == "backend." && k == :exec)) ||
            throw(ArgumentError("unknown solver option $(path)$(k); $(isempty(path) ? "sections" : "fields of $(chop(path))"): $(join(keys(ref), ", "))"))
    end
    function sub(k)
        r, u = get(ref, k, nothing), getfield(user, k)
        return r isa NamedTuple && u isa NamedTuple ? merge_config(r, u; path="$(path)$(k).") : u
    end
    return merge(ref, NamedTuple{keys(user)}(map(sub, keys(user))))
end

"""
    get_solver(; dim=2, kwargs...) -> ExplicitSolver or ImplicitSolver

The solver for a `dim`-dimensional problem: `get_default()` with the `kwargs` merged over it
(`merge_config`: per field, so a partial section keeps the other defaults; unknown top-level keys are
warned about, unknown nested fields are an error), `dtype.T0`/`bits` set from `dtype.precision` when
only the precision is given, the execution backend resolved and the kernel table (`cairn`) built.
"""
function get_solver(; dim::Number=2, kwargs...)
    ref = get_default()
    # Only keep keys in kwargs that are also in ref
    valids = intersect(keys(kwargs), keys(ref))
    user   = NamedTuple{Tuple(valids)}(getindex.(Ref(kwargs), valids))
    instr  = merge_config(ref, user)
    # a precision on its own brings its index and float types
    if haskey(user, :dtype) && user.dtype isa NamedTuple && haskey(user.dtype, :precision) && !haskey(user.dtype, :T0)
        o = get_option().dtype
        i = findfirst(==(user.dtype.precision), o.precision[2])
        isnothing(i) && throw(ArgumentError("dtype.precision = \"$(user.dtype.precision)\" is not one of $(o.precision[2])"))
        instr = merge(instr, (; dtype=merge(instr.dtype, (; T0=o.T0[2][i], bits=o.bits[2][i]))))
    end
    # Warn about unused keys
    unused = setdiff(keys(kwargs), keys(ref))
    if !isempty(unused)
        @warn join(vcat("miscellaneous kwargs:", "\n\t- ".*String.(unused)))
    end
    # Validate dimension
    if dim ∉ (1, 2, 3)
        error("Unsupported dimension: $dim")
    end
    # Verify plotting status from ui
    if !ElastoPlasm.self.ui.plot
        instr = merge(instr, (; plot = merge(instr.plot, (; status = false,)),))
    end
    return build_solver(instr, dim)
end

"""
    build_solver(instr, dim) -> ExplicitSolver or ImplicitSolver

The solver of a complete configuration `instr` (`get_default()`'s sections, merged): the execution backend selected
(`backend.select`, through UnifiedBackend, unless `backend.exec` is given) and the kernels built on it (`cairn`).
`get_solver` ends with it, and `load_simulation` rebuilds a saved solver with it from its stored configuration
(`solver_config`), so the device is picked on the machine that loads the file.
"""
function build_solver(instr::NamedTuple, dim::Integer)
    # Set execution backend
    if !haskey(instr[:backend],:exec)
        exec  = UnifiedBackend.select_execution_backend(UnifiedBackend.get_backend().exec, instr[:backend][:select]; distributed=instr[:backend][:distributed])
        length(exec) > 1 && @warn "$(length(exec)) devices selected: the solver runs on the first, $(first(values(exec)).name)"
        UnifiedBackend.device_wakeup!(first(values(exec)).handle)   # a GPU becomes the current device; nothing for a host
        instr = merge(instr, (; backend = merge(instr[:backend], (; exec = exec,)),))
    end
    # Add cairns to instr     
    cairn = (;
        ignite   = init_ignite(instr),
        mapsto   = init_mapsto(instr),
        update   = init_update(instr),
        implicit = init_implicit(instr),
    )
    instr = merge(instr, (; cairn = cairn,))
    # Infer the concrete solver type from instr[:solution] at construction time
    solver = if instr[:solution] == "explicit"
        ExplicitSolver
    elseif instr[:solution] == "implicit"
        ImplicitSolver
    else
        error("Unsupported solution: $(instr[:solution]) (expected \"explicit\" or \"implicit\")")
    end
    # Create Instruction type
    return solver{instr[:dtype][:T0]...,dim}(
        instr[:solution],
        instr[:dtype],
        instr[:basis],
        instr[:stab],
        instr[:bcs],
        instr[:grf],
        instr[:material],
        instr[:nonloc],
        instr[:plot],
        instr[:output],
        instr[:backend],
        instr[:cairn],
    )
end
