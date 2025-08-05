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

    T1, T2 = first(instr[:dtype].T0), last(instr[:dtype].T0)
    out = Mesh{T1,T2}(
        T1(mesh.dim), 
        T1.(mesh.nel), 
        T1.(mesh.nno), 
        T1(mesh.nn), 
        T2.(mesh.L), 
        T2.(mesh.h), 
        # nodal quantities
        T2.(mesh.x₀), 
        T2.(mesh.x), 
        T2.(mesh.mᵢ), 
        T2.(mesh.Mᵢⱼ),
        T2.(mesh.oobf), 
        T2.(mesh.D), 
        T2.(mesh.f), 
        T2.(mesh.a), 
        T2.(mesh.p), 
        T2.(mesh.v), 
        T2.(mesh.Δu), 
        T2.(mesh.ΔJ), 
        T2.(mesh.bij),
        # mesh-to-node topology
        T1.(mesh.e2n), 
        T1.(mesh.e2e), 
        T2.(mesh.xB), 
        # mesh boundary conditions
        T2.(mesh.bc)
    )
    return out
end