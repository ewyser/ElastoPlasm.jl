"""
    RFS(xp, zp, coh0, cohr, phi0, phir)

Generate random fields for cohesion and friction angle using exponential covariance.

# Arguments
- `xp`, `zp`: Coordinates of material points.
- `coh0`, `cohr`: Mean and residual cohesion values.
- `phi0`, `phir`: Mean and residual friction angle values.

# Returns
- `c`: Cohesion field (vector).
- `ϕ`: Friction angle field (vector).
"""
function RFS(xp,zp,coh0,cohr,phi0,phir)
    # parameters
    θx,θz     = 20.0,2.0
    β         = 45.0*π/180
    μc,σc     = coh0, coh0/5.0
    μϕ,σϕ     = phi0,phi0/10.0
    # vector format
    xp,zp,nmp = vec(xp),vec(zp),length(vec(xp))
    # relative distance
    Δx,Δz     = (xp.-xp'),zp.-zp'
    # exponential covariance matrix
    if β != 0.0
        C = real.(exp.(-sqrt.(complex.((( Δx.*cos(β).+Δz.*sin(β))./θx).^2+((-Δx.*cos(β).+Δz.*sin(β))./θz).^2))))
    else
        C = real.(exp.(-sqrt.(complex.((Δx./θx).^2+(Δz./θz).^2))))  
    end
    C[diagind(C)].= 1.0    
    cϕ   = cholesky(C).L*randn(Float64,nmp,2)
    p    = 0.5
    R    = [1.0 0.0;p sqrt(1.0-p^2)]
    cϕ   = R*cϕ'
    c    = μc.+σc.*cϕ[1,:]
    ϕ    = μϕ.+σϕ.*cϕ[2,:]

    p    = findall(x->x<=cohr, c)
    c[p].= cohr
	return c,ϕ
end

    #=
    #ρ  = exp.(-sqrt.(complex.((  Δx                     ./θx).^2+(  Δz                     ./θz).^2)))
    ρ  = exp.(-(complex.((Δx./θx).^2+(Δz./θz).^2)))
    C  = real.(ρ)
    Q  = eigvecs(C)
    Λ  = diagm(eigvals(C))
    c  = (Q*Λ.^(0.5)*randn(Float64,nmp)).+μ    
    =#
"""
    random_field(xp, σ, I, Nh, kₘ; covariance="gaussian") -> Vector

Zero-mean Gaussian random field of standard deviation `σ`, evaluated at the points `xp` (D×N; 2D
points are (x, z)), by the randomisation (spectral) method of Räss, Kolyukhin & Minakov (2019,
Computers & Geosciences 131:158-169), Eq. (3), with the phase taken at the points' own
coordinates (Algorithm 3). `I = (Ix, Iy, Iz)` are correlation lengths as defined by Eq. (6);
`Nh` is the number of harmonics.
- `covariance="gaussian"`: `C(r) = σ²exp(-r²/l_f²)`, `l_f = 2I/√π` (Eqs. 13-15, Algorithm 2);
  `kₘ` truncates the sampled dimensionless wave number. The paper's case is isotropic; unequal
  `I` scale each wave-vector component by its own `l_f` (an extension, equivalent to rescaling
  the axes, identical to Algorithm 2 when all `I` are equal).
- `covariance="exponential"`: `C(r) = σ²exp(-√Σ(rⱼ/Iⱼ)²)`, anisotropic as in the paper
  (Eqs. 5, 9-12, Algorithm 1); `kₘ` is unused.
"""
function random_field(xp::AbstractMatrix, σ, I, Nh, kₘ; covariance::String="gaussian")
    covariance ∈ ("gaussian","exponential") || error("covariance=\"$covariance\" is not implemented (\"gaussian\" or \"exponential\")")
    gauss = covariance == "gaussian"
    D     = size(xp,1)
    lf    = gauss ? 2.0 .* I ./ sqrt(π) : I     # per-axis length the wave vector is divided by
    f     = zeros(size(xp,2))
    for _ ∈ 1:Nh
        φ = 2.0*π*rand()
        k = 0.0
        if gauss                              # Alg. 2: p(k') ∝ k'² exp(-k'²/2) on [0,kₘ], k = k'√2
            while true
                k = kₘ*rand()
                rand()*2.0*exp(-1.0) < k*k*exp(-0.5*k*k) && break
            end
            k = sqrt(2.0)*k
        else                                  # Alg. 1: p(k) = 4k²/(π(1+k²)²)
            while true
                k = tan(0.5*π*rand())
                rand() < k*k/(1.0+k*k) && break
            end
        end
        θ   = acos(1.0-2.0*rand())
        V1  = k*sin(φ)*sin(θ)/lf[1]
        V2  = k*cos(φ)*sin(θ)/lf[2]
        V3  = k*cos(θ)/lf[3]
        a,b = randn(), randn()
        for p ∈ axes(xp,2)
            # phase at the point's own coordinates: x↔V1, (y↔V2 in 3D), vertical↔V3
            t     = xp[1,p]*V1 + xp[end,p]*V3 + (D == 3 ? xp[2,p]*V2 : 0.0)
            f[p] += a*sin(t) + b*cos(t)
        end
    end
    return (σ/sqrt(Nh)) .* f
end

"""
    get_cohesion(xp, mat, solver) -> Vector

Initial cohesion at the points `xp` (D×N): uniform `mat[:c0]`, or, with `solver.grf.status`,
`mat[:c0]` plus a Gaussian random field (`random_field`, covariance and parameters from `solver.grf`),
bounded below by the residual cohesion `mat[:cr]`.
"""
function get_cohesion(xp, mat, solver)
    solver.grf.status || return fill(mat[:c0], size(xp,2))
    g = solver.grf.param
    return max.(mat[:c0] .+ random_field(xp, g.σ, g.Iₓ, g.Nₕ, g.kₘ; covariance=solver.grf.covariance), mat[:cr])
end
