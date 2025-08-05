function setup_mps(mesh::Mesh{T1,T2},cmp::NamedTuple;define::Tuple=(nothing,nothing)) where {T1,T2}
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
    n0 = 0.0.*ones(nmp)
    l0 = ones(size(xp)).*0.5.*(mesh.h./ni)
    v0 = prod(2 .* l0; dims=1)
    ρ0 = ones(nmp) .* cmp[:ρ0]
    m  = (1.0 .- n0).*cmp[:ρ0].*v0
    # constructor
    mp = (
        ndim = mesh.dim,
        nmp  = nmp,
        vmax = zeros(mesh.dim),
        x    = copy(xp),
        z₀   = copy(xp[end,:]),
        n₀   = copy(n0),
        ℓ₀   = copy(l0), 
        ℓ    = copy(l0),
        Ω₀   = vec(copy(v0)),
        Ω    = vec(copy(v0)),
        s = (;
            u    = zeros(size(xp)), 
            v    = zeros(size(xp)),
            p    = zeros(size(xp)),

            ρ₀   = vec(copy(ρ0)),
            ρ    = vec(copy(ρ0)),
            m    = vec(copy(m)),
            c₀   = vec(copy(geom.coh0)),
            cᵣ   = vec(copy(geom.cohr)),
            ϕ    = vec(copy(geom.phi)),            
            Δλ   = zeros(nmp),
            ϵpII = zeros(2,nmp),
            ϵpV  = zeros(nmp), 
            ΔJ   = ones(nmp),
            J    = ones(nmp),
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
        ),
        l = (;

        ),
        # additional quantities
        ϕ∂ϕ  = zeros(mesh.nn,nmp ,mesh.dim+1   ),
        δnp  = zeros(mesh.nn,mesh.dim,nmp      ),
        # connectivity
        e2p  = spzeros(Int,nmp,mesh.nel[end]),
        p2p  = spzeros(Int,nmp,nmp),
        p2e  = zeros(Int,nmp),
        p2n  = zeros(Int,mesh.nn,nmp),
    )

    s = Solid{T1,T2}(
        T2.(mp.s.u)    ,
        T2.(mp.s.v)    ,
        T2.(mp.s.p)    ,
        # mechanical properties
        T2.(mp.s.ρ₀)   ,
        T2.(mp.s.ρ)    ,
        T2.(mp.s.m)    ,
        T2.(mp.s.c₀)   ,
        T2.(mp.s.cᵣ)   ,
        T2.(mp.s.ϕ)    ,
        T2.(mp.s.Δλ)   ,
        T2.(mp.s.ϵpII) ,
        T2.(mp.s.ϵpV)  ,
        T2.(mp.s.ΔJ)   ,
        T2.(mp.s.J)    ,
        # tensor in voigt notation
        T2.(mp.s.σᵢ)   ,
        T2.(mp.s.τᵢ)   ,
        # tensor in matrix notation
        T2.(mp.s.δᵢⱼ)  ,
        T2.(mp.s.∇vᵢⱼ) ,
        T2.(mp.s.∇uᵢⱼ) ,
        T2.(mp.s.ΔFᵢⱼ) ,
        T2.(mp.s.Fᵢⱼ)  ,
        T2.(mp.s.Bᵢⱼ)  ,
        T2.(mp.s.ϵᵢⱼ)  ,
        T2.(mp.s.ωᵢⱼ)  ,
        T2.(mp.s.σJᵢⱼ) ,
    )
    l = Liquid{T1,T2}(

    )
    out = Point{T1,T2}(
        # general information
        T1(mp.ndim) ,
        T1(mp.nmp)  ,
        T2.(mp.vmax) ,
        # basis-related quantities
        T2.(mp.ϕ∂ϕ)  ,
        T2.(mp.δnp)  ,
        # connectivity
        T1.(mp.e2p)  ,
        T1.(mp.p2p)  ,
        T1.(mp.p2e)  ,
        T1.(mp.p2n)  ,
        # material point properties
        T2.(mp.x)    ,
        T2.(mp.ℓ₀)   ,
        T2.(mp.ℓ)    ,
        T2.(mp.Ω₀)   ,
        T2.(mp.Ω)    ,
        # solid phase
        s       ,
        # liquid phase
        l       ,
    )
    return out 
end
