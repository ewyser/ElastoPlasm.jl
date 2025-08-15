"""
    init_shpfun(dim::Number, instr::NamedTuple; what::String="nothing")

Initialize shape function and topology kernels for the MPM algorithm.

# Arguments
- `dim::Number`: Spatial dimension (1, 2, or 3).
- `instr::NamedTuple`: Instruction/configuration dictionary.
- `what::String`: (Optional) Additional selector for kernel type.

# Returns
- Named tuple of kernel functions for topology, shape function, and delta function.
"""
function init_shpfun(dim::Number,instr::NamedTuple;what::String="nothing")
    kernel3,kernel4 = nothing,nothing
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
        kernel3 = Δnp_nd(CPU())   
        kernel4 = nothing  
    elseif instr[:fwrk][:trsfr] == "apic"
        kernel3 = nothing
        kernel4 = Dij_nd(CPU()) 
    end
    return (;tplgy! = kernel1, ϕ∂ϕ! = kernel2, Δₙₚ! = kernel3, Dᵢⱼ! = kernel4,)
end
"""
    shpfun(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, instr::NamedTuple) where {T1,T2}

Initialize and compute shape functions and topological relations for material points.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function shpfun(mpts::Point{T1,T2},mesh::Mesh{T1,T2},instr::NamedTuple) where {T1,T2} 
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    instr[:cairn][:shpfun].tplgy!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    # initialize shapefunctions
    mpts.ϕ∂ϕ .= T2(0.0)
    # calculate shape functions
    instr[:cairn][:shpfun].ϕ∂ϕ!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    # calculate identity shape functions
    if instr[:fwrk][:trsfr] == "tpic"
        instr[:cairn][:shpfun].Δₙₚ!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    elseif instr[:fwrk][:trsfr] == "apic"
        instr[:cairn][:shpfun].Dᵢⱼ!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    end
    return nothing
end