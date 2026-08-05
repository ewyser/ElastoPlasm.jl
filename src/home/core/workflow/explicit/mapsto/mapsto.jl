"""
    init_mapsto(dim::Number, instr::NamedTuple)

Initialize mapping, solving and transfering kernels for MPM cycle based on dimension, material phase and instruction set.

# Arguments
- `dim::Number`: Spatial dimension (1, 2, or 3).
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `Dict`: Dictionary of mapping and augmentation kernels.
"""
function init_mapsto(instr::NamedTuple; mapsto::Dict = Dict(:map => Dict{Symbol,Cairn}(),:augm => Dict{Symbol,Cairn}()))
    if instr.fwrk.deform == "finite"
        mapsto[:map][:σᵢ!] = transform(CPU())
    end
    if instr.fwrk.trsfr == "std"
        mapsto[:map][:p2n!] = std_p2n(CPU())
    elseif instr.fwrk.trsfr == "tpic"
        mapsto[:map][:p2n!] = tpic_p2n(CPU())
    elseif instr.fwrk.trsfr == "apic"
        mapsto[:map][:p2n!] = apic_p2n(CPU())
        mapsto[:map][:Bᵢⱼ!] = Bij(CPU())
    else
        return throw(ArgumentError("$(instr.fwrk.trsfr) is an unsupported transfer scheme"))
    end
    mapsto[:map][:solve!] = euler(CPU())
    mapsto[:map][:n2p!]   = picflip_n2p(CPU())
    if instr.fwrk.musl 
        mapsto[:augm][:p2n!]   = augm_p2n(CPU())
        mapsto[:augm][:solve!] = augm_solve(CPU())
        return (; map = (; mapsto[:map]...), augm = (; mapsto[:augm]...))
    else
        return (; map = (; mapsto[:map]...), augm = nothing)
    end
end
"""
    mapsto(mpts::Point{T1,T2,D,E,R},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2},dt::T2,instr::NamedTuple) where {T1,T2,D,E,R}

Resolution of mechanical problem: project material points to nodes, solve, and map back.

# Arguments
- `mpts::Point{T1,T2,D,E,R}`: Material point data structure.
- `mesh::MeshSolidPhase{T1,T2,D}`: Mesh data structure for solid phase.
- `g::Vector{T2}`: Gravity vector.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function mapsto(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},g::Vector{T2},dt::T2,instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R}
    # get cauchy stress 
    if instr.fwrk.deform == "finite"
        instr.cairn.mapsto.map.σᵢ!(ndrange=mpts.nmp,mpts);sync(CPU())
    end
    # reset nodal quantities
    fill!(mesh.s.m  ,T2(0.0))
    fill!(mesh.s.mv  ,T2(0.0))
    fill!(mesh.s.oobf,T2(0.0))
    fill!(mesh.s.a   ,T2(0.0))
    fill!(mesh.s.v   ,T2(0.0))
    # mapping to mesh
    instr.cairn.mapsto.map.p2n!(mpts,mesh,g; ndrange=mpts.nmp);sync(CPU())
    # solve Eulerian momentum equation
    instr.cairn.mapsto.map.solve!(mesh,dt,T2(instr.fwrk.damping); ndrange=mesh.prprt.nno[end]);sync(CPU())
    # maps back solution to material point
    instr.cairn.mapsto.map.n2p!(mpts,mesh,dt,T2(instr.fwrk.C_pf); ndrange=mpts.nmp);sync(CPU())
    # (if musl) reproject nodal velocities
    if instr.fwrk.musl
        # reset nodal quantities
        fill!(mesh.s.mv,T2(0.0))
        fill!(mesh.s.v ,T2(0.0))
        # accumulate material point contributions
        instr.cairn.mapsto.augm.p2n!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # solve for nodal incremental displacement
        instr.cairn.mapsto.augm.solve!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
    end
    # (for APIC) compute Bᵢⱼ for material points
    if instr.fwrk.trsfr == "apic"
        instr.cairn.mapsto.map.Bᵢⱼ!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
    end
    return nothing
end
