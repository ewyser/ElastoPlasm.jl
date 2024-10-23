function pointSetup(meD,cmParam,instr;define::Tuple=(nothing,nothing))
    # non-dimensional constant                                                   
    if meD.nD == 2 
        nstr = 3 
    elseif meD.nD == 3 
        nstr = 6 
    end
    # material geometry
    ni,nmp,geom = define
    # scalars & vectors
    if meD.nD == 2
        l0 = ones(typeD,nmp,meD.nD).*0.5.*(meD.h[1]./ni)  
        v0 = ones(typeD,nmp).*(2.0.*l0[:,1].*2.0.*l0[:,2])
    elseif meD.nD == 3
        l0 = ones(typeD,nmp,meD.nD).*0.5.*(meD.h[1]./ni)                
        v0 = ones(typeD,nmp).*(2.0.*l0[:,1].*2.0.*l0[:,2].*2.0.*l0[:,3])
    end
    m = cmParam[:ρ0].*v0
    # constructor
    mpD = (
        nD   = meD.nD,
        nmp  = nmp,
        x    = copy(geom.xp),
        u    = zeros(typeD,nmp,meD.nD), 
        v    = zeros(typeD,nmp,meD.nD),
        p    = zeros(typeD,nmp,meD.nD),
        ℓ₀   = copy(l0), 
        ℓ    = copy(l0),
        Ω₀   = copy(v0),
        Ω    = copy(v0),
        m    = copy(m),
        c₀   = copy(geom.coh0),
        cᵣ   = copy(geom.cohr),
        ϕ    = copy(geom.phi),            
        Δλ   = zeros(typeD,nmp  ),
        ϵpII = zeros(typeD,nmp,2),
        ϵpV  = zeros(typeD,nmp), 
        ΔJ   = ones(typeD,nmp),
        J    = ones(typeD,nmp),
        # plot quantity
        z₀   = copy(geom.xp[:,end]),
        # tensor in matrix notation
        I    = Matrix(1.0I,meD.nD,meD.nD    ),
        ∇vᵢⱼ = zeros(typeD,meD.nD,meD.nD,nmp),
        ∇uᵢⱼ = zeros(typeD,meD.nD,meD.nD,nmp),
        ΔFᵢⱼ = zeros(typeD,meD.nD,meD.nD,nmp),
        Fᵢⱼ  = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        ϵᵢⱼ  = zeros(typeD,meD.nD,meD.nD,nmp),
        ωᵢⱼ  = zeros(typeD,meD.nD,meD.nD,nmp),
        σJᵢⱼ = zeros(typeD,meD.nD,meD.nD,nmp),
        Bᵢⱼ  = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        # tensor in voigt notation
        σᵢ   = zeros(typeD,nstr,nmp),
        τᵢ   = zeros(typeD,nstr,nmp),
        # additional quantities
        ϕ∂ϕ  = zeros(typeD,meD.nn,nmp ,meD.nD+1   ),
        δnp  = zeros(typeD,meD.nn,meD.nD,nmp      ),
        # connectivity
        e2p  = spzeros(Int64,nmp,meD.nel[end]),
        p2p  = spzeros(Int64,nmp,nmp),
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,meD.nn,nmp),
    )
    return mpD 
end
