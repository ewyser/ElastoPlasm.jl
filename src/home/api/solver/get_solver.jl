export get_solver

function get_solver(; dim::Number=2, kwargs...)
    ref = get_default()
    # Only keep keys in kwargs that are also in ref
    valids = intersect(keys(kwargs), keys(ref))
    user   = NamedTuple{Tuple(valids)}(getindex.(Ref(kwargs), valids))
    instr  = merge(ref, user)
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
    # 
    if instr[:perf][:status]
        instr = merge(instr, (; strain = merge(instr[:strain], (; deform = "infinitesimal",)),))
        instr = merge(instr, (; nonloc = merge(instr[:nonloc], (; status = false,)),))
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
        instr[:strain],
        instr[:stab],
        instr[:bcs],
        instr[:grf],
        instr[:plast],
        instr[:nonloc],
        instr[:plot],
        instr[:perf],
        instr[:backend],
        instr[:cairn],
    )
end
