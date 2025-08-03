function setup_mesh(nel,L,instr)
    # geometry                                               
    if instr[:basis][:which] == "bsmpm"
        ndim,nn,h = length(L),4^length(L),min.(L./nel,L./4)
    else
        ndim,nn,h = length(L),4^length(L),L./nel #L,h,ndim,nn = getinfo(L,nel)
    end

    if instr[:basis][:ghost]
        #@info "Init Eulerian mesh & adding ghosts"
        buffer = 2.0.*h
    else
        #@info "Init Eulerian mesh"
        buffer = 0.0.*h            
    end
    # mesh & boundary conditions
    xn,nel,nno = get_coords(ndim,nn,L,h;ghosts=buffer)
    bc,xB      = get_bc(xn,h,nno,ndim;ghosts=buffer)
    # constructor 
    mesh = (;
        dim  = ndim,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h, # mᵢ Mᵢⱼ
        # nodal quantities
        x₀   = vec(minimum(xn,dims=2)     ),
        x    = xn                          ,
        mᵢ   = zeros(nno[end]             ), # lumped mass vector
        Mᵢⱼ  = zeros(nno[end],nno[end]    ), # consistent mass matrix
        oobf = zeros(ndim,nno[end]        ),
        D    = zeros(ndim,nno[end]        ),
        f    = zeros(ndim,nno[end]        ),
        a    = zeros(ndim,nno[end]        ),
        p    = zeros(ndim,nno[end]        ),
        v    = zeros(ndim,nno[end]        ),
        Δu   = zeros(ndim,nno[end]        ),
        ΔJ   = zeros(ndim,nno[end]        ),
        bij  = zeros(ndim,ndim,nno[end]   ),
        # mesh-to-node topology
        e2n  = e2n(ndim,nno,nel,nn        ),
        e2e  = e2e(ndim,nno,nel,nn,h,instr),
        xB   = xB                          ,
        # mesh boundary conditions
        bc   = bc                          ,
    )
    #=
    println("")
    println("Material Point Setup:")
    for key in keys(mesh)
        println("$(key): $(typeof(mesh[key]))")
    end
    =#
    T1, T2     = Int64              , Float64
    A3, A5, A7 = AbstractArray{T1,1}, AbstractArray{T1,2}, AbstractArray{T1,3}
    A2, A4, A6 = AbstractArray{T2,1}, AbstractArray{T2,2}, AbstractArray{T2,3}
    out = Mesh{T1,T2,A3,A5,A7,A2,A4,A6}(
        mesh.dim, 
        mesh.nel, 
        mesh.nno, 
        mesh.nn, 
        mesh.L, 
        mesh.h, 
        # nodal quantities
        mesh.x₀, 
        mesh.x, 
        mesh.mᵢ, 
        mesh.Mᵢⱼ,
        mesh.oobf, 
        mesh.D, 
        mesh.f, 
        mesh.a, 
        mesh.p, 
        mesh.v, 
        mesh.Δu, 
        mesh.ΔJ, 
        mesh.bij,
        # mesh-to-node topology
        mesh.e2n, 
        mesh.e2e, 
        mesh.xB, 
        # mesh boundary conditions
        mesh.bc
    )
    return out
end