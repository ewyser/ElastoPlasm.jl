function setup_mps(mesh,cmp;define::Tuple=(nothing,nothing))
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
    l0 = ones(size(xp)).*0.5.*(mesh.h./ni)
    v0 = prod(2 .* l0; dims=1)
    m  = cmp[:ρ0].*v0
    # constructor
    mp = (
        ndim = mesh.dim,
        nmp  = nmp,
        vmax = zeros(mesh.dim),
        x    = copy(xp),
        u    = zeros(size(xp)), 
        v    = zeros(size(xp)),
        p    = zeros(size(xp)),
        ℓ₀   = copy(l0), 
        ℓ    = copy(l0),
        Ω₀   = vec(copy(v0)),
        Ω    = vec(copy(v0)),
        m    = vec(copy(m)),
        c₀   = vec(copy(geom.coh0)),
        cᵣ   = vec(copy(geom.cohr)),
        ϕ    = vec(copy(geom.phi)),            
        Δλ   = zeros(nmp),
        ϵpII = zeros(2,nmp),
        ϵpV  = zeros(nmp), 
        ΔJ   = ones(nmp),
        J    = ones(nmp),
        # plot quantity
        z₀   = copy(xp[end,:]),
        # tensor in matrix notation
        δᵢⱼ  = Matrix(1.0I,mesh.dim,mesh.dim), 
        ∇vᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
        ∇uᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
        ΔFᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
        Fᵢⱼ  = repeat(Matrix(1.0I,mesh.dim,mesh.dim),1,1,nmp),
        Bᵢⱼ  = repeat(Matrix(1.0I,mesh.dim,mesh.dim),1,1,nmp),
        ϵᵢⱼ  = zeros(mesh.dim,mesh.dim,nmp),
        ωᵢⱼ  = zeros(mesh.dim,mesh.dim,nmp),
        σJᵢⱼ = zeros(mesh.dim,mesh.dim,nmp),
        # tensor in voigt notation
        σᵢ   = zeros(nstr,nmp),
        τᵢ   = zeros(nstr,nmp),
        # additional quantities
        ϕ∂ϕ  = zeros(mesh.nn,nmp ,mesh.dim+1   ),
        δnp  = zeros(mesh.nn,mesh.dim,nmp      ),
        # connectivity
        e2p  = spzeros(Int,nmp,mesh.nel[end]),
        p2p  = spzeros(Int,nmp,nmp),
        p2e  = zeros(Int,nmp),
        p2n  = zeros(Int,mesh.nn,nmp),
    )
    #=
    println("")
    println("Material Point Setup:")
    for key in keys(mp)
        println("$(key): $(typeof(mp[key]))")
    end
    =#
    T1, T2     = Int64              , Float64
    A3, A5, A7 = AbstractArray{T1,1}, AbstractArray{T1,2}, AbstractArray{T1,3}
    A2, A4, A6 = AbstractArray{T2,1}, AbstractArray{T2,2}, AbstractArray{T2,3}
    out = Point{T1,T2,A3,A5,A7,A2,A4,A6}(
        mp.ndim ,
        mp.nmp  ,
        mp.vmax ,
        mp.x    ,
        mp.u    ,
        mp.v    ,
        mp.p    ,
        mp.ℓ₀   ,
        mp.ℓ    ,
        mp.Ω₀   ,
        mp.Ω    ,
        mp.m    ,
        mp.c₀   ,
        mp.cᵣ   ,
        mp.ϕ    ,
        mp.Δλ   ,
        mp.ϵpII ,
        mp.ϵpV  ,
        mp.ΔJ   ,
        mp.J    ,
        # plot quantity
        mp.z₀   ,
        # tensor in matrix notation
        mp.δᵢⱼ  ,
        mp.∇vᵢⱼ ,
        mp.∇uᵢⱼ ,
        mp.ΔFᵢⱼ ,
        mp.Fᵢⱼ  ,
        mp.Bᵢⱼ  ,
        mp.ϵᵢⱼ  ,
        mp.ωᵢⱼ  ,
        mp.σJᵢⱼ ,
        # tensor in voigt notation
        mp.σᵢ   ,
        mp.τᵢ   ,
        # additional quantities
        mp.ϕ∂ϕ  ,
        mp.δnp  ,
        # connectivity
        mp.e2p  ,
        mp.p2p  ,
        mp.p2e  ,
        mp.p2n  ,
    )
    return mp 
end
