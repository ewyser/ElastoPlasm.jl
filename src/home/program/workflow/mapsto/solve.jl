"""
    euler(mesh::MeshSolidPhase{T1,T2},dt::T2,η::T2) where {T2}

Solve the Eulerian momentum equation for solid phase on the mesh with viscous damping.

# Arguments
- `mesh::MeshSolidPhase{T1,T2}`: Mesh data structure for solid phase.
- `dt::T2`: Time step.
- `η`: Damping coefficient.

# Returns
- Updates mesh fields in-place.
"""
@views @kernel inbounds = true function euler(mesh::MeshSolidPhase{T1,T2},dt::T2,η::T2) where {T1,T2}
    no = @index(Global)
    if no ≤ mesh.prprt.nno[end]
        if iszero(mesh.mᵢ[no])
            nothing         
        else
            for dim ∈ 1:mesh.prprt.dim
                # apply boundary contidions
                if mesh.bcs.status[dim,no]
                    mesh.a[dim,no] = T2(0.0)                                         
                    mesh.v[dim,no] = T2(0.0)                                         
                else
                    # cache mass node & norm of out-of-balance force
                    mᵢ = (T2(1.0)/mesh.mᵢ[no])
                    # calculate damping
                    D  = η*norm(mesh.oobf[:,no])*sign(mesh.mv[dim,no]*mᵢ)       #(2,)
                    if (abs(mesh.mv[dim,no]*mᵢ)) ≥ T2(1.0e-3)
                        mesh.oobf[dim,no] = mesh.oobf[dim,no]-D                 #(2,)
                    end
                    # forward euler solution
                    mesh.a[dim,no] = mesh.oobf[dim,no]*mᵢ                        #(2,)
                    mesh.v[dim,no] = (mesh.mv[dim,no]+dt*mesh.oobf[dim,no])*mᵢ #(2,)  
                end
            end
        end
    end
end
"""
    euler(mesh::MeshThermalPhase{T1,T2},dt::T2,η::T2) where {T2}

Solve the Eulerian momentum equation for thermal phase on the mesh with viscous damping.

# Arguments
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure for thermal phase.
- `dt::T2`: Time step.
- `η`: Damping coefficient.

# Returns
- Updates mesh fields in-place.
"""
@views @kernel inbounds = true function euler(mesh::MeshThermalPhase{T1,T2},dt::T2) where {T1,T2}
    no = @index(Global)
    if no≤mesh.prprt.nno[end] 
        if iszero(mesh.cᵢ[no])
            nothing         
        else
            # apply boundary contidions
            if mesh.bcs.status[1,no]
                mesh.T[no] = T2(20.0)
            else
                # cache mass node & norm of out-of-balance force
                cᵢ = (T2(1.0)/mesh.cᵢ[no])
                # forward euler solution
                mesh.dT[no] = mesh.oobq[no]*cᵢ                        #(1,)
                mesh.T[no]  = (mesh.mcT[no]+dt*(mesh.oobq[no]))*cᵢ    #(1,)
            end
        end
    end
end
"""
    solve(mesh::MeshSolidPhase{T1,T2}, dt::T2, instr::NamedTuple)

Solve the mesh momentum equation for solid phase using the backend-agnostic kernel.

# Arguments
- `mesh::MeshSolidPhase{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates mesh fields in-place.
"""
function solve(mesh::MeshSolidPhase{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # viscous damping
    η = T2(instr.fwrk.damping)
    # initialize
    fill!(mesh.a,T2(0.0))
    fill!(mesh.v,T2(0.0))
    # solve momentum equation on the mesh using backend-agnostic kernel
    instr[:cairn][:mapsto][:map].solve!(mesh,dt,η; ndrange=mesh.prprt.nno[end]);sync(CPU())
    return nothing
end
"""
    solve(mesh::MeshThermalPhase{T1,T2}, dt::T2, instr::NamedTuple)

Solve the mesh momentum equation for thermal phase using the backend-agnostic kernel.

# Arguments
- `mesh::MeshThermalPhase{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates mesh fields in-place.
"""
function solve(mesh::MeshThermalPhase{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # initialize
    fill!(mesh.dT,T2(0.0))
    fill!(mesh.T,T2(0.0))
    # solve momentum equation on the mesh using backend-agnostic kernel
    instr[:cairn][:mapsto][:map].solve!(mesh,dt; ndrange=mesh.prprt.nno[end]);sync(CPU())
    return nothing
end