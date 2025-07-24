@views @kernel inbounds = true function euler(mesh,dt,η)
    no = @index(Global)
    for dim ∈ 1:mesh.dim
        if no≤mesh.nno[end] 
            if iszero(mesh.mn[no])
                m = 0.0                                                          #(1,)                       
            else
                m = (1.0/mesh.mn[no])*mesh.bc[dim,no]                            #(1,)
            end
            mesh.Dn[dim,no] = η*norm(mesh.oobf[:,no])*sign(mesh.pn[dim,no]*m)    #(2,)
            mesh.fn[dim,no] = mesh.oobf[dim,no]-mesh.Dn[dim,no]                  #(2,)
            mesh.an[dim,no] = mesh.fn[dim,no]*m                                  #(2,)
            mesh.vn[dim,no] = (mesh.pn[dim,no]+dt*mesh.fn[dim,no])*m             #(2,)   
        end
    end
end
@views function solve(mesh,dt,instr)
    # viscous damping
    η      = 0.1
    # initialize
    mesh.fn.= 0.0
    mesh.an.= 0.0
    mesh.vn.= 0.0
    # solve momentum equation on the mesh using backend-agnostic kernel
    instr[:cairn][:mapsto][:map].solve!(ndrange=mesh.nno[end],mesh,dt,η);sync(CPU())
    return nothing
end