@views @kernel inbounds = true function euler(mesh,dt,η)
    no = @index(Global)
    for dim ∈ 1:mesh.dim
        if no≤mesh.nno[end] 
            if iszero(mesh.mᵢ[no])
                m = 0.0                                                          #(1,)                       
            else
                m = (1.0/mesh.mᵢ[no])*mesh.bc[dim,no]                            #(1,)
            end
            mesh.D[dim,no] = η*norm(mesh.oobf[:,no])*sign(mesh.p[dim,no]*m)    #(2,)
            mesh.f[dim,no] = mesh.oobf[dim,no]-mesh.D[dim,no]                  #(2,)
            mesh.a[dim,no] = mesh.f[dim,no]*m                                  #(2,)
            mesh.v[dim,no] = (mesh.p[dim,no]+dt*mesh.f[dim,no])*m             #(2,)   
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