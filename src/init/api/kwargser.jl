function kwargser(type::Symbol,kwargs::Any;dim::Number = 2)
    ref  = require(type)
    key0 = collect(keys(kwargs))
    
    k,v  = [],[]
    if type == :instr
        for (i,key) ∈ enumerate(collect(keys(ref)))
            push!(k,key)
            if haskey(kwargs,key)
                push!(v,kwargs[key])
                filter!(x->x≠key,key0)
            else
                push!(v,ref[key])
            end
        end
        # display @warn message
        if !isempty(key0)
            @warn join(vcat("miscellaneous kwargs for require():","\n\t- ".*String.(key0)))
        end
        # zip & set precision
        instr = Dict(zip(k,v))
        bits  = pop!(instr,:bits); delete!(instr,:bits)
        if bits == 64
            instr[:dtype] = (;T0=(Int64,Float64),bits=64,precision="double")
        elseif bits == 32
            instr[:dtype] = (;T0=(Int32,Float32),bits=32,precision="single")
        end
        # add cairns (abstract kernels) to instr set
        instr[:cairn] = (;
            tplgy!  = init_shpfun(dim,instr[:basis];what="tplgy!"),
            ϕ∂ϕ!    = init_shpfun(dim,instr[:basis];what="ϕ∂ϕ!"  ),
            p2n!    = init_mapsto(dim,instr[:trsfr];what="p2n!"  ),
            n2p!    = init_mapsto(dim,instr[:trsfr];what="n2p!"  ),
            augm!   = init_double(instr[:trsfr]),
            solveEuler! = init_solve(),
            deform! = init_deformation(instr),
            ΔJn!    = init_volumetric(;what="ΔJn!"),
            ΔJs!    = init_volumetric(;what="ΔJs!"),
            ΔJp!    = init_volumetric(;what="ΔJp!"),
            update! = init_domain(dim,instr[:basis]),
        )
        return instr
    else
        for (i,var) ∈ enumerate(key)
            if haskey(kwargs,var)
                push!(values,kwargs[var])
                if var == :select
                    selection = kwargs[var]
                end
            else
                push!(values,ref[var])
                if var == :select
                    selection = ref[var]
                end
            end
        end
        return Dict(zip(key,value)) #return (;zip(Tuple(Symbol(x) for x in key),value)...)
    end
end