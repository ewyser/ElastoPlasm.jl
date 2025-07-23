function pointSetup(mesh,cmParam,instr;define::Tuple=(nothing,nothing))
    # non-dimensional constant                                                   
    if mesh.dim == 2 
        nstr = 3 
    elseif mesh.dim == 3 
        nstr = 6 
    end
    # material geometry
    ni,nmp,geom = define
    xp = geom.xp 
    # scalars & vectors
    if mesh.dim == 2
        l0 = ones(size(xp)).*0.5.*(mesh.h[1]./ni)  
        v0 = ones(nmp).*(2.0.*l0[1,:].*2.0.*l0[2,:])
    elseif mesh.dim == 3
        l0 = ones(size(xp)).*0.5.*(mesh.h[1]./ni)                
        v0 = ones(nmp).*(2.0.*l0[1,:].*2.0.*l0[2,:].*2.0.*l0[3,:])
    end
    m = cmParam[:ρ0].*v0
    # constructor
    mp = (
        ndim   = mesh.dim,
        nmp  = nmp,
        x    = copy(xp),
        u    = zeros(size(xp)), 
        v    = zeros(size(xp)),
        p    = zeros(size(xp)),
        ℓ₀   = copy(l0), 
        ℓ    = copy(l0),
        Ω₀   = copy(v0),
        Ω    = copy(v0),
        m    = copy(m),
        c₀   = copy(geom.coh0),
        cᵣ   = copy(geom.cohr),
        ϕ    = copy(geom.phi),            
        Δλ   = zeros(nmp),
        ϵpII = zeros(size(xp)),
        ϵpV  = zeros(nmp), 
        ΔJ   = ones(nmp),
        J    = ones(nmp),
        # plot quantity
        z₀   = copy(xp[end,:]),
        # tensor in matrix notation
        δᵢⱼ  = Matrix(1.0I,mesh.dim,mesh.dim    ), 
        ∇vᵢⱼ = zeros(typeD,mesh.dim,mesh.dim,nmp),
        ∇uᵢⱼ = zeros(typeD,mesh.dim,mesh.dim,nmp),
        ΔFᵢⱼ = zeros(typeD,mesh.dim,mesh.dim,nmp),
        Fᵢⱼ  = repeat(Matrix(1.0I,mesh.dim,mesh.dim),1,1,nmp),
        Bᵢⱼ  = repeat(Matrix(1.0I,mesh.dim,mesh.dim),1,1,nmp),
        ϵᵢⱼ  = zeros(typeD,mesh.dim,mesh.dim,nmp),
        ωᵢⱼ  = zeros(typeD,mesh.dim,mesh.dim,nmp),
        σJᵢⱼ = zeros(typeD,mesh.dim,mesh.dim,nmp),
        # tensor in voigt notation
        σᵢ   = zeros(typeD,nstr,nmp),
        τᵢ   = zeros(typeD,nstr,nmp),
        # additional quantities
        ϕ∂ϕ  = zeros(typeD,mesh.nn,nmp ,mesh.dim+1   ),
        δnp  = zeros(typeD,mesh.nn,mesh.dim,nmp      ),
        # connectivity
        e2p  = spzeros(Int64,nmp,mesh.nel[end]),
        p2p  = spzeros(Int64,nmp,nmp),
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,mesh.nn,nmp),
    )
    return mp 
end
