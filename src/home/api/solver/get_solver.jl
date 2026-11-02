
"""
    exec_backend(x) -> KernelAbstractions backend

The backend a solver runs on, for the solver or its configuration (`instr`): its platform's (`backend.exec`, the
UnifiedBackend `Platform` `build_solver` selected), `CPU()` for the host. The kernels are built on it (`init_ignite`,
`init_mapsto`, `init_update`, `init_implicit`), the workflows synchronize on it and `to_device` moves a problem to it:
the solver is where the backend is read, never the arrays.
"""
function exec_backend(x)
    return x.backend.exec.backend
end

"""
    solver_config(solver) -> NamedTuple

A solver's configuration, the sections `build_solver` takes, with its run settings: what `save_simulation` keeps of a
solver, without the selected platform (`backend.exec`) or the kernels (`cairn`), which `load_simulation` rebuilds where
the file is loaded. `stored_solver(solver)` is the solver holding only it (`cairn = nothing`).
"""
function solver_config(s::Solver)
    return (; solution=solution_name(s), dtype=s.dtype, basis=s.basis, stab=s.stab, nonloc=s.nonloc, plot=s.plot,
              output=s.output, backend=(; select=s.backend.select), run=s.run)
end
function stored_solver(s::Solver{T1,T2,D}) where {T1,T2,D}
    c = solver_config(s)
    return Solver{T1,T2,D}(s.solution, c.dtype, c.basis, c.stab, c.nonloc, c.plot, c.output, c.backend, c.run, nothing)
end

"""
    merge_config(ref, user; path="") -> NamedTuple

`user`'s settings over `ref`'s, section by section: nested sections are merged field by field
(`(; basis=(; which="gimpm"))` changes only `basis.which`), anything else (values, vectors such as
`plot.what`) replaces the default. A field `ref` does not have is an error naming the valid ones
(`backend.exec`, set by `build_solver`, is accepted).
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
    Solver(; dim=2, sections…) -> Solver

A solver for `dim`-dimensional cases, written in code (`Solver(case)` reads one from a case's settings):
`get_default()` with `sections` merged over it (`solver_instr`), the execution backend resolved and the kernel table
(`cairn`) built; the run settings are `RUN`'s. The physics (`material`, `grf`, `bcs`) is a case's, so it is refused
here (`ArgumentError`).
"""
function Solver(; dim::Integer=2, kwargs...)
    refuse_physics(kwargs)
    return build_solver(solver_instr(; dim, kwargs...), dim)
end

"""
    solver_instr(; dim=2, kwargs...) -> NamedTuple

A complete solver configuration: `get_default()` with the `kwargs` merged over it (`merge_config`: per field, so a
partial section keeps the other defaults; unknown top-level keys are warned about, unknown nested fields are an
error), `dtype.T0`/`bits` set from `dtype.precision` when only the precision is given, and plotting off when the
session has no plots.
"""
function solver_instr(; dim::Number=2, kwargs...)
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
    return instr
end

"""
    build_solver(instr, dim) -> Solver

The solver of a complete configuration `instr` (`get_default()`'s numerical sections, merged; with its run settings as
`run`, else `RUN`'s, `solver_run`): the execution backend selected (`backend.select`, through UnifiedBackend, unless
`backend.exec` is given) and the kernels built on it (`cairn`). Both `Solver` constructors end with it, and
`load_simulation` rebuilds a saved solver with it from its stored configuration (`solver_config`), so the device is
picked on the machine that loads the file.
"""
function build_solver(instr::NamedTuple, dim::Integer)
    # Set execution backend
    if !haskey(instr[:backend],:exec)
        exec = UnifiedBackend.select_platform(instr[:backend][:select])
        UnifiedBackend.activate!(exec)   # a GPU becomes the current device
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
    return Solver{instr[:dtype][:T0]...,dim}(
        get_solution(instr[:solution]),                     # the solution's marker, what the workflows dispatch on
        instr[:dtype],
        instr[:basis],
        instr[:stab],
        instr[:nonloc],
        instr[:plot],
        instr[:output],
        instr[:backend],
        haskey(instr, :run) ? instr[:run] : solver_run((;)),
        instr[:cairn],
    )
end
