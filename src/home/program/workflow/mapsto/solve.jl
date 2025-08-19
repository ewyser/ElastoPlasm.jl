"""
    euler(mesh, dt::T2, η) where {T2}

Solve the Eulerian momentum equation on the mesh with viscous damping.

# Arguments
- `mesh`: Mesh data structure.
- `dt::T2`: Time step.
- `η`: Damping coefficient.

# Returns
- Updates mesh fields in-place.
"""
@views @kernel inbounds = true function euler(::Val{:Solid},mesh::Mesh{T1,T2},dt::T2,η::T2) where {T1,T2}
    no = @index(Global)
    if no≤mesh.nno[end]
        if iszero(mesh.s.mᵢ[no])
            nothing         
        else
            for dim ∈ 1:mesh.dim
                # apply boundary contidions
                if mesh.bcs.status[dim,no]
                    mesh.s.a[dim,no] = T2(0.0)                                         
                    mesh.s.v[dim,no] = T2(0.0)                                         
                else
                    # cache mass node & norm of out-of-balance force
                    mᵢ = (T2(1.0)/mesh.s.mᵢ[no])
                    # calculate damping
                    D  = η*norm(mesh.s.oobf[:,no])*sign(mesh.s.mv[dim,no]*mᵢ)       #(2,)
                    if (abs(mesh.s.mv[dim,no]*mᵢ)) ≥ T2(1.0e-3)
                        mesh.s.oobf[dim,no] = mesh.s.oobf[dim,no]-D                 #(2,)
                    end
                    # forward euler solution
                    mesh.s.a[dim,no] = mesh.s.oobf[dim,no]*mᵢ                        #(2,)
                    mesh.s.v[dim,no] = (mesh.s.mv[dim,no]+dt*mesh.s.oobf[dim,no])*mᵢ #(2,)  
                end
            end
        end
    end
end
@views @kernel inbounds = true function euler(::Val{:Thermal},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    no = @index(Global)
    if no≤mesh.nno[end] 
        if iszero(mesh.t.cᵢ[no])
            nothing         
        else
            # apply boundary contidions
            if mesh.bcs.status[dim,no]
                #mesh.T[no] = T2(0.0)                                         
            else
                # cache mass node & norm of out-of-balance force
                cᵢ = (T2(1.0)/mesh.t.cᵢ[no])
                # forward euler solution
                mesh.t.T[no]+= dt*(mesh.t.oobq[dim,no]*cᵢ)                        #(2,)
            end
        end
    end
end
"""
    solve(mesh, dt::T2, instr::NamedTuple) where {T2}

Solve the mesh momentum equation using the backend-agnostic kernel.

# Arguments
- `mesh`: Mesh data structure.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates mesh fields in-place.
"""
@views function solve(mesh::Mesh{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # viscous damping
    η      = T2(instr.fwrk.damping)
    # initialize
    fill!(mesh.s.a,T2(0.0))
    fill!(mesh.s.v,T2(0.0))
    # solve momentum equation on the mesh using backend-agnostic kernel
    instr[:cairn][:mapsto][:map].solve!(Val(:Solid),mesh,dt,η; ndrange=mesh.nno[end]);sync(CPU())
    return nothing
end