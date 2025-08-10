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
@views @kernel inbounds = true function euler(mesh,dt::T2,η) where {T2}
    no = @index(Global)
    if no≤mesh.nno[end] 
        if iszero(mesh.mᵢ[no])
            nothing         
        else
            # cache mass node & norm of out-of-balance force
            mᵢ,oobf = (T2(1.0)/mesh.mᵢ[no]),norm(mesh.oobf[:,no])
            for dim ∈ 1:mesh.dim
                # apply boundary contidions
                if mesh.bcs.status[dim,no]
                    mesh.a[dim,no] = T2(0.0)                                         
                    mesh.v[dim,no] = T2(0.0)                                         
                else
                    # calculate damping
                    mesh.D[dim,no] = η*oobf*sign(mesh.p[dim,no]*mᵢ)                  #(2,)
                    # forward euler solution
                    mesh.f[dim,no] = mesh.oobf[dim,no]-mesh.D[dim,no]                #(2,)
                    mesh.a[dim,no] = mesh.f[dim,no]*mᵢ                               #(2,)
                    mesh.v[dim,no] = (mesh.p[dim,no]+dt*mesh.f[dim,no])*mᵢ           #(2,)  
                end
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
@views function solve(mesh,dt::T2,instr::NamedTuple) where {T2}
    # viscous damping
    η      = T2(instr.fwrk.damping)
    # initialize
    mesh.f.= T2(0.0)
    mesh.a.= T2(0.0)
    mesh.v.= T2(0.0)
    # solve momentum equation on the mesh using backend-agnostic kernel
    instr[:cairn][:mapsto][:map].solve!(ndrange=mesh.nno[end],mesh,dt,η);sync(CPU())
    return nothing
end