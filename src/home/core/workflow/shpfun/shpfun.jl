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
function init_shpfun(instr::NamedTuple; shpfun::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}(),)
    # topology function
    shpfun[:tplgy!] = p2e2n(CPU())
    # shape function
    if instr[:basis][:which] == "bsmpm"
        shpfun[:ϕ∂ϕ!] = bsmpm(CPU())
    elseif instr[:basis][:which] == "gimpm"
        shpfun[:ϕ∂ϕ!] = gimpm(CPU())
    elseif instr[:basis][:which] == "smpm"
        shpfun[:ϕ∂ϕ!] = smpm(CPU())   
    else
        return throw(ArgumentError("$(instr[:basis][:which]) is not a supported shape function basis"))
    end
    if instr[:fwrk][:trsfr] == "tpic"
        shpfun[:Δₙₚ!] = Δnp_nd(CPU()) 
    elseif instr[:fwrk][:trsfr] == "apic"
        shpfun[:Dᵢⱼ!] = Dij_nd(CPU()) 
    end    
    return (;shpfun...)
end
"""
    shpfun(mpts::Point{T1,T2,E,R}, mesh::Mesh{T1,T2,D}, instr::NamedTuple) where {T1,T2,E,R,D}

Initialize and compute shape functions and topological relations for material points.

# Arguments
- `mpts::Point{T1,T2,E,R}`: Material point data structure.
- `mesh::Mesh{T1,T2,D}`: Mesh data structure.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function shpfun(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R} 
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    instr.cairn.shpfun.tplgy!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    # initialize shapefunctions
    fill!(mpts.ϕ∂ϕ,T2(0.0))
    # calculate shape functions
    instr.cairn.shpfun.ϕ∂ϕ!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    # calculate identity shape functions
    if instr.fwrk.trsfr == "tpic"
        instr.cairn.shpfun.Δₙₚ!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    elseif instr.fwrk.trsfr == "apic"
        instr.cairn.shpfun.Dᵢⱼ!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    end
    return nothing
end