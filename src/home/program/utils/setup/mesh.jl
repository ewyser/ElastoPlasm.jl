function meshSetup(nel,L,instr)
    # geometry                                               
    L,h,ndim,nn = getinfo(L,nel)
    if instr[:basis][:ghost]
        @info "Init Eulerian mesh & adding ghosts"
        buffer = 2.0.*h
    else
        @info "Init Eulerian mesh"
        buffer = 0.0.*h            
    end
    # mesh 
    xn,nel,nno = get_coords(ndim,nn,L,h;ghosts=buffer)
    # boundary conditions
    bc,xB      = get_bc(xn,h,nno,ndim;ghosts=buffer)
    # constructor 
    mesh = (
        dim  = ndim,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h,
        # nodal quantities
        x₀   = vec(minimum(xn,dims=2)),
        xn   = xn,
        mn   = zeros(nno[end]              ), # lumped mass vector
        Mn   = zeros(nno[end],nno[end]     ), # consistent mass matrix
        oobf = zeros(ndim,nno[end]         ),
        Dn   = zeros(ndim,nno[end]         ),
        fn   = zeros(ndim,nno[end]         ),
        an   = zeros(ndim,nno[end]         ),
        pn   = zeros(ndim,nno[end]         ),
        vn   = zeros(ndim,nno[end]         ),
        Δun  = zeros(ndim,nno[end]         ),
        ΔJn  = zeros(ndim,nno[end]         ),
        bijn = zeros(ndim,ndim,nno[end]    ),
        # mesh-to-node topology
        e2n  = e2n(ndim,nno,nel,nn),
        e2e  = e2e(ndim,nno,nel,nn,h,instr),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return mesh
end