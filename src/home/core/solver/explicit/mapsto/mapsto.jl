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
    if instr.strain.deform == "finite"
        mapsto[:map][:σᵢ!] = transform(CPU())
    end
    if instr.transfer.trsfr == "std"
        mapsto[:map][:p2n!] = std_p2n(CPU())
    elseif instr.transfer.trsfr == "tpic"
        mapsto[:map][:p2n!] = tpic_p2n(CPU())
    elseif instr.transfer.trsfr == "apic"
        mapsto[:map][:p2n!] = apic_p2n(CPU())
        mapsto[:map][:Bᵢⱼ!] = Bij(CPU())
    else
        return throw(ArgumentError("$(instr.transfer.trsfr) is an unsupported transfer scheme"))
    end
    mapsto[:map][:solve!] = euler(CPU())
    mapsto[:map][:n2p!]   = picflip_n2p(CPU())
    if instr.transfer.musl 
        mapsto[:augm][:p2n!]   = augm_p2n(CPU())
        mapsto[:augm][:solve!] = augm_solve(CPU())
        return (; map = (; mapsto[:map]...), augm = (; mapsto[:augm]...))
    else
        return (; map = (; mapsto[:map]...), augm = nothing)
    end
end
"""
    mapsto(mpts::Point{T1,T2,D,E},mesh::MeshSolidPhase{T1,T2,D},g::Vector{T2},dt::T2,instr::NamedTuple) where {T1,T2,D,E}

Resolution of mechanical problem: project material points to nodes, solve, and map back.

# Arguments
- `mpts::Point{T1,T2,D,E}`: Material point data structure.
- `mesh::MeshSolidPhase{T1,T2,D}`: Mesh data structure for solid phase.
- `g::Vector{T2}`: Gravity vector.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function mapsto(mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},g::Vector{T2},dt::T2,solver::ExplicitSolver{T1,T2}) where {T1,T2}
    # get cauchy stress
    if solver.strain.deform == "finite"
        solver.cairn.mapsto.map.σᵢ!(ndrange=mpts.nmp,mpts);sync(CPU())
    end
    # reset nodal quantities
    fill!(mesh.s.m  ,T2(0.0))
    fill!(mesh.s.mv  ,T2(0.0))
    fill!(mesh.s.oobf,T2(0.0))
    fill!(mesh.s.a   , zero(eltype(mesh.s.a)))
    fill!(mesh.s.v   , zero(eltype(mesh.s.v)))
    # mapping to mesh
    solver.cairn.mapsto.map.p2n!(mpts,mesh,basis,g; ndrange=mpts.nmp);sync(CPU())
    # solve Eulerian momentum equation
    solver.cairn.mapsto.map.solve!(mesh,dt,T2(solver.stab.damping); ndrange=mesh.prprt.nno[end]);sync(CPU())
    # maps back solution to material point
    solver.cairn.mapsto.map.n2p!(mpts,mesh,basis,dt,T2(solver.transfer.C_pf); ndrange=mpts.nmp);sync(CPU())
    # (if musl) reproject nodal velocities
    if solver.transfer.musl
        # reset nodal quantities
        fill!(mesh.s.mv, T2(0.0))
        fill!(mesh.s.v , zero(eltype(mesh.s.v)))
        # accumulate material point contributions
        solver.cairn.mapsto.augm.p2n!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
        # solve for nodal incremental displacement
        solver.cairn.mapsto.augm.solve!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
    end
    # (for APIC) compute Bᵢⱼ for material points
    if solver.transfer.trsfr == "apic"
        solver.cairn.mapsto.map.Bᵢⱼ!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
    end
    return nothing
end

"""
    mapsto(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},basis::Basis{T1,T2},dt::T2,solver::ExplicitSolver{T1,T2}) where {T1,T2}

Resolution of the thermal (heat-conduction) problem: project material-point heat
capacity/temperature to nodes, solve the nodal heat-balance equation, and map the
updated temperature back to material points. Mirrors the solid-phase `mapsto` above,
minus the gravity/finite-strain/APIC concerns that don't apply to the thermal phase
(only the `std` transfer scheme has a thermal `p2n!` kernel — see `p2n/std.jl`).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure for the thermal phase.
- `basis::Basis{T1,T2}`: Topology container linking `mesh`/`mpts`.
- `dt::T2`: Time step.
- `solver::ExplicitSolver{T1,T2}`: Solver instance.

# Returns
- `nothing`. Updates fields in-place.
"""
function mapsto(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},basis::Basis{T1,T2},dt::T2,solver::ExplicitSolver{T1,T2}) where {T1,T2}
    # reset nodal quantities
    fill!(mesh.cᵢ  ,T2(0.0))
    fill!(mesh.mcT ,T2(0.0))
    fill!(mesh.oobq,T2(0.0))
    fill!(mesh.dT  ,T2(0.0))
    # mapping to mesh
    solver.cairn.mapsto.map.p2n!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
    # solve nodal heat-balance equation
    solver.cairn.mapsto.map.solve!(mesh,dt; ndrange=mesh.prprt.nno[end]);sync(CPU())
    # maps back solution to material point
    solver.cairn.mapsto.map.n2p!(mpts,mesh,basis,dt,T2(solver.transfer.C_pf); ndrange=mpts.nmp);sync(CPU())
    # (if musl) reproject nodal temperature
    if solver.transfer.musl
        fill!(mesh.mcT,T2(0.0))
        solver.cairn.mapsto.augm.p2n!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
        solver.cairn.mapsto.augm.solve!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
    end
    return nothing
end
