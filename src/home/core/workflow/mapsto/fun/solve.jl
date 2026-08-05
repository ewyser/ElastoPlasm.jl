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
@views @kernel inbounds = true function euler(mesh::Mesh{T1,T2},dt::T2,η::T2) where {T1,T2}
    no = @index(Global)
    if no ≤ mesh.prprt.nno[end]
        if iszero(mesh.s.m[no])
            nothing         
        else
            m⁻¹ = (T2(1.0)/mesh.s.m[no])
            for dim ∈ 1:mesh.prprt.dim
                # apply boundary contidions
                if mesh.s.bcs.status[dim,no]
                    mesh.s.a[dim,no] = T2(0.0)                                         
                    mesh.s.v[dim,no] = T2(0.0)                                         
                else
                    # calculate damping
                    D  = η*norm(mesh.s.oobf[:,no])*sign(mesh.s.mv[dim,no]*m⁻¹)      
                    if (abs(mesh.s.mv[dim,no]*m⁻¹)) ≥ T2(1.0e-3)
                        mesh.s.oobf[dim,no] = mesh.s.oobf[dim,no]-D                 
                    end
                    # forward euler solution
                    mesh.s.a[dim,no] = mesh.s.oobf[dim,no]*m⁻¹                       
                    mesh.s.v[dim,no] = (mesh.s.mv[dim,no]+dt*mesh.s.oobf[dim,no])*m⁻¹  
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