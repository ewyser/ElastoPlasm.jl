"""
    init_mapsto(dim::Number, instr::NamedTuple)

Initialize mapping and transfer kernels for the MPM algorithm based on dimension and instruction set.

# Arguments
- `dim::Number`: Spatial dimension (1, 2, or 3).
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `Dict`: Dictionary of mapping and augmentation kernels.
"""
function init_mapsto(dim::Number,instr::NamedTuple) 
    mapsto = Dict(:map => Dict(),)
    if instr[:fwrk][:deform] == "finite"
        mapsto[:map][:σᵢ!] = transform(CPU())
    end
    if instr[:fwrk][:trsfr] == "std"
        if dim == 1
            mapsto[:map][:p2n!] = std_1d_p2n(CPU()) 
        elseif dim == 2
            mapsto[:map][:p2n!] = std_2d_p2n(CPU())
        elseif dim == 3
            mapsto[:map][:p2n!] = std_3d_p2n(CPU())
        end
    elseif instr[:fwrk][:trsfr] == "tpic"
        if dim == 1
            mapsto[:map][:p2n!] = tpic_1d_p2n(CPU())
        elseif dim == 2
            mapsto[:map][:p2n!] = tpic_2d_p2n(CPU())
        elseif dim == 3
            mapsto[:map][:p2n!] = tpic_3d_p2n(CPU())
        end
    elseif instr[:fwrk][:trsfr] == "apic"
        if dim == 1
            nothing # APIC is not yet defined for 1D
        elseif dim == 2
            mapsto[:map][:p2n!] = apic_2d_p2n(CPU())
        elseif dim == 3
            mapsto[:map][:p2n!] = apic_3d_p2n(CPU())
        end
        mapsto[:map][:Bᵢⱼ!] = Bij(CPU())
    else
        return throw(ArgumentError("$(instr[:fwrk][:trsfr]) is an unsupported transfer scheme"))
    end
    mapsto[:map][:solve!] = euler(CPU())
    mapsto[:map][:n2p!]   = picflip_n2p(CPU())
    if instr[:fwrk][:musl] 
        mapsto[:augm]          = Dict()
        mapsto[:augm][:p2n!]   = augm_momentum(CPU())
        mapsto[:augm][:solve!] = augm_velocity(CPU())
    end
    return Dict(:map => (; mapsto[:map]...), :augm => (; mapsto[:augm]...))
end
"""
    mapsto(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, g::Vector{T2}, dt::T2, instr::NamedTuple) where {T1,T2}

Perform a full MPM update: project material points to nodes, solve, and map back.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `g::Vector{T2}`: Gravity vector.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function mapsto(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2},g::Vector{T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # maps material point to node
    p2n(mpts,mesh,g,instr)
    # solve Eulerian momentum equation
    solve(mesh,dt,instr)
    # maps back solution to material point
    n2p(mpts,mesh,dt,instr)
    return nothing
end