@views @kernel inbounds = true function euler(mesh,dt,η)
    no = @index(Global)
    if no≤mesh.nno[end] 
        if iszero(mesh.mᵢ[no])
            nothing         
        else
            # cache mass node & norm of out-of-balance force
            mᵢ,oobf = (1.0/mesh.mᵢ[no]),norm(mesh.oobf[:,no])
            for dim ∈ 1:mesh.dim
                # calculate damping
                mesh.D[dim,no] = η*oobf*sign(mesh.p[dim,no]*mᵢ)                  #(2,)
                # forward euler solution
                mesh.f[dim,no] = mesh.oobf[dim,no]-mesh.D[dim,no]                #(2,)
                mesh.a[dim,no] = mesh.f[dim,no]*mᵢ                               #(2,)
                mesh.v[dim,no] = (mesh.p[dim,no]+dt*mesh.f[dim,no])*mᵢ           #(2,)  
                # apply boundary contidions
                if mesh.bc[dim,no] == 0.0
                    mesh.a[dim,no] = 0.0                                         
                    mesh.v[dim,no] = 0.0                                         
                end
            end
        end
    end
end
@views function solve(mesh,dt,instr)
    # viscous damping
    η      = 0.1
    # initialize
    mesh.f.= 0.0
    mesh.a.= 0.0
    mesh.v.= 0.0
    # solve momentum equation on the mesh using backend-agnostic kernel
    instr[:cairn][:mapsto][:map].solve!(ndrange=mesh.nno[end],mesh,dt,η);sync(CPU())
    return nothing
end