"""
    mapsto(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2},g::Vector{T2},dt::T2,instr::NamedTuple)

Resolution of mechanical problem: project material points to nodes, solve, and map back.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshSolidPhase{T1,T2}`: Mesh data structure for solid phase.
- `g::Vector{T2}`: Gravity vector.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function mapsto(mpts::Point{T1,T2},mesh::MeshSolidPhase{T1,T2},g::Vector{T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # get cauchy stress 
    if instr[:fwrk][:deform] == "finite"
        instr[:cairn][:mapsto][:map].σᵢ!(ndrange=mpts.nmp,mpts);sync(CPU())
    end
    # reset nodal quantities
    fill!(mesh.mᵢ  ,T2(0.0))
    fill!(mesh.mv  ,T2(0.0))
    fill!(mesh.oobf,T2(0.0))
    fill!(mesh.a   ,T2(0.0))
    fill!(mesh.v   ,T2(0.0))
    # mapping to mesh
    instr[:cairn][:mapsto][:map].p2n!(mpts,mesh,g; ndrange=mpts.nmp);sync(CPU())
    # solve Eulerian momentum equation
    instr[:cairn][:mapsto][:map].solve!(mesh,dt,T2(instr.fwrk.damping); ndrange=mesh.prprt.nno[end]);sync(CPU())
    # maps back solution to material point
    instr[:cairn][:mapsto][:map].n2p!(mpts,mesh,dt,T2(instr[:fwrk][:C_pf]); ndrange=mpts.nmp);sync(CPU())
    # (if musl) reproject nodal velocities
    if instr[:fwrk][:musl]
        # reset nodal quantities
        fill!(mesh.mv,T2(0.0))
        fill!(mesh.v ,T2(0.0))
        # accumulate material point contributions
        instr[:cairn][:mapsto][:augm].p2n!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # solve for nodal incremental displacement
        instr[:cairn][:mapsto][:augm].solve!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
    end
    # (for APIC) compute Bᵢⱼ for material points
    if instr[:fwrk][:trsfr] == "apic"
        instr[:cairn][:mapsto][:map].Bᵢⱼ!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
    end
    return nothing
end
"""
    mapsto(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},dt::T2,instr::NamedTuple)

Resolution of thermal problem: project material points to nodes, solve, and map back.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure for thermal phase.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function mapsto(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # reset nodal quantities
    fill!(mesh.cᵢ  ,T2(0.0))
    fill!(mesh.mcT ,T2(0.0))
    fill!(mesh.oobq,T2(0.0))
    fill!(mesh.dT  ,T2(0.0))
    fill!(mesh.T   ,T2(0.0))
    # mapping to mesh
    instr[:cairn][:mapsto][:map].p2n!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
    # solve Eulerian thermal equation
    instr[:cairn][:mapsto][:map].solve!(mesh,dt; ndrange=mesh.prprt.nno[end]);sync(CPU())
    # maps back solution to material point
    instr[:cairn][:mapsto][:map].n2p!(mpts,mesh,dt,T2(instr[:fwrk][:C_pf]); ndrange=mpts.nmp);sync(CPU())
    # (if musl) reproject nodal velocities
    if instr[:fwrk][:musl]
        # reset nodal quantities
        fill!(mesh.mcT,T2(0.0))
        fill!(mesh.T  ,T2(0.0))
        # accumulate material point contributions
        instr[:cairn][:mapsto][:augm].p2n!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # solve for nodal temperature
        instr[:cairn][:mapsto][:augm].solve!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
    end
    return nothing
end
"""
    mapsto(mpts::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2},dt::T2,instr::NamedTuple)

Resolution of thermo-poro-mechanical problem: project material points to nodes, solve, and map back.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `g::Vector{T2}`: Gravity vector.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function mapsto(mpts::Point{T1,T2},mesh::Mesh{T1,T2},g::Vector{T2},dt::T2,instr::NamedTuple) where {T1,T2}

    return nothing
end
