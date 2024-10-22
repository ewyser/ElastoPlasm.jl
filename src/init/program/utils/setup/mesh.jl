function meshSetup(nel,L,instr;ghost::Bool=false)
    @info "init Eulerian mesh geometry"
    # geometry                                               
    L,h,nD,nn    = getinfo(L,nel)
    buffer       = 0.0.*h
    if instr[:basis] == "gimpm" && ghost
        buffer = 2.0.*h
    end
    # mesh 
    xn,tn,nel,nno = getcoords(nD,nn,L,h;ghosts=buffer)
    # boundary conditions
    bc,xB         = getbc(xn,h,nno,nD;ghosts=buffer)
    # constructor 
    meD = (
        nD   = nD,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h,
        # nodal quantities
        x₀   = minimum(xn,dims=1),
        xn   = xn,
        tn   = Int64.(tn),
        mn   = zeros(nno[end]            ), # lumped mass vector
        Mn   = zeros(nno[end],nno[end]   ), # consistent mass matrix
        oobf = zeros(nno[end],nD         ),
        Dn   = zeros(nno[end],nD         ),
        fn   = zeros(nno[end],nD         ),
        an   = zeros(nno[end],nD         ),
        pn   = zeros(nno[end],nD         ),
        vn   = zeros(nno[end],nD         ),
        Δun  = zeros(nno[end],nD         ),
        ΔJn  = zeros(nno[end],nD         ),
        bijn = zeros(nD      ,nD,nno[end]),
        # mesh-to-node topology
        e2n  = e2n(nD,nno,nel,nn),
        e2e  = e2e(nD,nno,nel,nn,h,instr),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return meD
end