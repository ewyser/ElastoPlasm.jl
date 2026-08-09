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
function init_ignite(instr::NamedTuple; ignite::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}(),)
    # topology function
    ignite[:tplgy!] = p2e2n(CPU())
    if instr[:fwrk][:trsfr] == "apic"
        ignite[:Dᵢⱼ!] = Dij_nd(CPU()) 
    end    
    return (;ignite...)
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
function ignite(mpts::Point{T1,T2,D,<:AbstractBasis,E,R},mesh::Mesh{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R} 
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    instr.cairn.ignite.tplgy!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    # calculate identity shape functions
    if instr.fwrk.trsfr == "apic"
        instr.cairn.ignite.Dᵢⱼ!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    end
    return nothing
end