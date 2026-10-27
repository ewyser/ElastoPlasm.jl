export get_solver

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
    # Set execution backend
    if !haskey(instr[:backend],:exec)
        exec  = select_execution_backend(ElastoPlasm.self.bckd, instr[:backend][:select]; distributed=instr[:backend][:distributed], prompt=ElastoPlasm.self.bckd.prompt)
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
