function init_shpfun(dim::Number,basis::NamedTuple;what::String="nothing")
    # topology function
    if dim == 1
        tplgy! = p2e2n1D!(CPU())
    elseif dim == 2
        tplgy! = p2e2n2D!(CPU())
    elseif dim == 3
        tplgy! = p2e2n3D!(CPU())
    end
    # shape function
    if basis[:which] == "bsmpm"
        if dim == 1
            ϕ∂ϕ! = bsmpm1D(CPU())    
        elseif dim == 2
            ϕ∂ϕ! = bsmpm2D(CPU())
        elseif dim == 3
            ϕ∂ϕ! = bsmpm3D(CPU())
        end
    elseif basis[:which] == "gimpm"
        if dim == 1
            ϕ∂ϕ! = gimpm1D(CPU())    
        elseif dim == 2
            ϕ∂ϕ! = gimpm2D(CPU())
        elseif dim == 3
            ϕ∂ϕ! = gimpm3D(CPU())
        end
    elseif basis[:which] == "smpm"
        if dim == 1
            ϕ∂ϕ! = smpm1D(CPU())    
        elseif dim == 2
            ϕ∂ϕ! = smpm2D(CPU())
        elseif dim == 3
            ϕ∂ϕ! = smpm3D(CPU())
        end
    else
        return throw(ArgumentError("$(basis[:which]) is not a supported shape function basis"))
    end
    return (;tplgy! = tplgy!, ϕ∂ϕ! = ϕ∂ϕ!,)
end
function shpfun!(mpD,meD,instr)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    instr[:cairn][:shpfun!].tplgy!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    # initialize shapefunctions
    mpD.ϕ∂ϕ .= 0.0
    # calculate shape functions
    instr[:cairn][:shpfun!].ϕ∂ϕ!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end