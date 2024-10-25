function init_shpfun(dim::Number,basis::NamedTuple;what::String="nothing")
    # topology function
    if dim == 1
        kernel1 = p2e2n1D!(CPU())
    elseif dim == 2
        kernel1 = p2e2n2D!(CPU())
    elseif dim == 3
        kernel1 = p2e2n3D!(CPU())
    end
    # shape function
    if basis[:which] == "bsmpm"
        if dim == 1
            kernel2 = bsmpm1D(CPU())    
        elseif dim == 2
            kernel2 = bsmpm2D(CPU())
        elseif dim == 3
            kernel2 = bsmpm3D(CPU())
        end
    elseif basis[:which] == "gimpm"
        if dim == 1
            kernel2 = gimpm1D(CPU())    
        elseif dim == 2
            kernel2 = gimpm2D(CPU())
        elseif dim == 3
            kernel2 = gimpm3D(CPU())
        end
    elseif basis[:which] == "smpm"
        if dim == 1
            kernel2 = smpm1D(CPU())    
        elseif dim == 2
            kernel2 = smpm2D(CPU())
        elseif dim == 3
            kernel2 = smpm3D(CPU())
        end
    else
        return throw(ArgumentError("$(basis[:which]) is not a supported shape function basis"))
    end
    return (;tplgy! = kernel1, ϕ∂ϕ! = kernel2,)
end
function shpfun(mpD,meD,instr)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    instr[:cairn][:shpfun].tplgy!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    # initialize shapefunctions
    mpD.ϕ∂ϕ .= 0.0
    # calculate shape functions
    instr[:cairn][:shpfun].ϕ∂ϕ!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end