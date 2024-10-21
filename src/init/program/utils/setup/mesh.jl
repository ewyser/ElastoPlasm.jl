function meshSetup(nel,L,instr;ghost::Bool=false)
    @info "init Eulerian mesh geometry"
    # geometry                                               
    L,h,nD       = meshGeom(L,nel)
    if instr[:basis] == "gimpm" && ghost
        buffer = 2.0.*h
    else
        buffer = 0.0.*h
    end
    # mesh 
    x,t,nn,nel,nno = meshCoord(nD,L,h;ghosts=buffer)
    # boundary conditions
    bc,xB        = meshBCs(x,h,nno,nD;ghosts=buffer)
    if nD>1
        minC = minimum(x,dims=1)
    else
        minC = minimum(x)
    end
    # constructor 
    meD = (
        nD   = nD,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h,
        # nodal quantities
        x₀   = minC,
        xn   = x,
        tn   = Int64.(t),
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
        bn   = zeros(nD      ,nD,nno[end]),
        # mesh-to-node topology
        e2n  = e2n(nD,nno,nel,nn),
        e2e  = e2e(nD,nno,nel,nn,h,instr),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return meD
end