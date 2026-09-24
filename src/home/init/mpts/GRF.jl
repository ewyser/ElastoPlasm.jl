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
    grf_gauss(xp, σ, I, Nh, kₘ) -> Vector

Zero-mean Gaussian random field with Gaussian covariance `C(r) = σ²exp(-r²/l_f²)`,
`l_f = 2I/√π`, evaluated at the points `xp` (D×N; 2D points are (x, z)), by the randomisation
(spectral) method of Räss, Kolyukhin & Minakov (2019, Computers & Geosciences 131:158-169),
Eqs. (3), (13)-(15) and Algorithms 2-3, with the phase taken at the points' own coordinates.
`I = (Ix, Iy, Iz)` are per-axis correlation lengths: equal values reproduce the paper's isotropic
Algorithm 2 exactly; unequal values scale each wave-vector component by its own `l_f` (an
extension of the paper, equivalent to rescaling the axes). `Nh` harmonics, `kₘ` maximum sampled
(dimensionless) wave number.
"""
function grf_gauss(xp::AbstractMatrix, σ, I, Nh, kₘ)
    D  = size(xp,1)
    lf = 2.0 .* I ./ sqrt(π)
    f  = zeros(size(xp,2))
    for _ ∈ 1:Nh
        φ = 2.0*π*rand()
        k = 0.0
        while true                        # rejection sampling of p(k') ∝ k'² exp(-k'²/2)
            k = kₘ*rand()
            rand()*2.0*exp(-1.0) < k*k*exp(-0.5*k*k) && break
        end
        θ   = acos(1.0-2.0*rand())
        V1  = sqrt(2.0)*k*sin(φ)*sin(θ)/lf[1]
        V2  = sqrt(2.0)*k*cos(φ)*sin(θ)/lf[2]
        V3  = sqrt(2.0)*k*cos(θ)/lf[3]
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
`mat[:c0]` plus a Gaussian random field (`grf_gauss`, parameters from `solver.grf.param`),
bounded below by the residual cohesion `mat[:cr]`.
"""
function get_cohesion(xp, mat, solver)
    solver.grf.status || return fill(mat[:c0], size(xp,2))
    solver.grf.covariance == "gaussian" || error("grf.covariance=\"$(solver.grf.covariance)\" is not implemented (only \"gaussian\")")
    g = solver.grf.param
    return max.(mat[:c0] .+ grf_gauss(xp, g.σ, g.Iₓ, g.Nₕ, g.kₘ), mat[:cr])
end

function GRFS_exp(xl,coh0,cohr,ni,Δx)
    # =====================================================================
    # GRFS: Gaussian Random Field Simulator - exponential covariance
    # 
    # Copyright (C) 2019  Ludovic Raess, Dmitriy Kolyukhin and Alexander Minakov.
    # 
    # GRFS is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.
    # 
    # GRFS is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    # 
    # You should have received a copy of the GNU General Public License
    # along with GRFS. If not, see <http://www.gnu.org/licenses/>.
    # =====================================================================
    # physics
    sf  = 5e3                        # standard deviation
    If  = [2.5,2.5,2.5]              # correlation lengths in [x,y,z]
    # numerics
    Nh  = 5000                       # inner parameter, number of harmonics
    nx  = size(xl,2)                 # numerical grid resolution in x
    ny  = size(xl,3)                 # numerical grid resolution in y
    nz  = size(xl,1)                 # numerical grid resolution in z
    dx  = 0.5*Δx/ni            # numerical grid step size in x
    dy  = 0.5*Δx/ni            # numerical grid step size in y
    dz  = 0.5*Δx/ni            # numerical grid step size in z
    # preprocessing
    C   = sf/sqrt(Nh)
    coh = zeros(Float64,nz,nx,ny)
    tmp = zeros(Float64,nz,nx,ny)
    # action
    for ih in 1:Nh
        r  = rand(1)
        fi = 2.0*π*r[1]

        # Gaussian spectrum
        flag = true;
    
        r = rand(1)
        r = r[1]
        k = tan(pi*0.5*r)
        while flag
            r = rand(1)
            r = r[1]
            k = tan(pi*0.5*r)
            d = (k*k)/(1.0+(k*k))
            r = rand(1)
            r = r[1]
            if r<d
                flag = false
            end
        end 
    
    
        r = rand(1)
        r = r[1]
        theta = acos(1-2*r)
        V1 = k*sin(fi)*sin(theta)/If[1]
        V2 = k*cos(fi)*sin(theta)/If[2]
        V3 = k*cos(theta)/If[3]
        r  = randn(1)
        a  = r[1]
        r  = randn(1)
        b  = r[1]
    
        for iz in 1:nz
            for iy in 1:ny
                for ix in 1:nx
                tmp[iz,ix,iy] = dx*(ix-0.5)*V1 + dy*(iy-0.5)*V2 + dz*(iz-0.5)*V3;
                end
            end
        end 
        coh = coh .+ a.*sin.(tmp) .+ b.*cos.(tmp);
    end
    coh .= coh0.+C.*coh;
    p   = findall(x->x<=cohr, coh)
    coh[p] .= cohr
    return coh
end