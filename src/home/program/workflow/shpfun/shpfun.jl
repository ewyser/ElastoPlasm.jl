function init_shpfun(dim::Number,instr::Dict;what::String="nothing")
    # topology function
    if dim == 1
        kernel1 = p2e2n_1d(CPU())
    elseif dim == 2
        kernel1 = p2e2n_2d(CPU())
    elseif dim == 3
        kernel1 = p2e2n_3d(CPU())
    end
    # shape function
    if instr[:basis][:which] == "bsmpm"
        if dim == 1
            kernel2 = bsmpm_1d(CPU())    
        elseif dim == 2
            kernel2 = bsmpm_2d(CPU())
        elseif dim == 3
            kernel2 = bsmpm_3d(CPU())
        end
    elseif instr[:basis][:which] == "gimpm"
        if dim == 1
            kernel2 = gimpm_1d(CPU())    
        elseif dim == 2
            kernel2 = gimpm_2d(CPU())
        elseif dim == 3
            kernel2 = gimpm_3d(CPU())
        end
    elseif instr[:basis][:which] == "smpm"
        if dim == 1
            kernel2 = smpm_1d(CPU())    
        elseif dim == 2
            kernel2 = smpm_2d(CPU())
        elseif dim == 3
            kernel2 = smpm_3d(CPU())
        end
    else
        return throw(ArgumentError("$(instr[:basis][:which]) is not a supported shape function basis"))
    end
    if instr[:fwrk][:trsfr] == "tpic"
        if dim == 1
            kernel3 = δ_1d(CPU())    
        elseif dim == 2
            kernel3 = δ_2d(CPU())
        elseif dim == 3
            kernel3 = δ_3d(CPU())
        end
    else
        kernel3 = nothing
    end
    return (;tplgy! = kernel1, ϕ∂ϕ! = kernel2, δ! = kernel3)
end
function shpfun(mp,mesh,instr)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    instr[:cairn][:shpfun].tplgy!(mp,mesh; ndrange=(mp.nmp));sync(CPU())
    # initialize shapefunctions
    mp.ϕ∂ϕ .= 0.0
    # calculate shape functions
    instr[:cairn][:shpfun].ϕ∂ϕ!(mp,mesh; ndrange=(mp.nmp));sync(CPU())
    # calculate identity shape functions
    if instr[:fwrk][:trsfr] == "tpic"
        instr[:cairn][:shpfun].δ!(mp,mesh; ndrange=(mp.nmp));sync(CPU())
    end
    return nothing
end