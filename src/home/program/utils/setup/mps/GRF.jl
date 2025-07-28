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
function GRFS_gauss(xl,coh0,cohr,ni,Δx)
    # =====================================================================
    # GRFS: Gaussian Random Field Simulator - Gaussian covariance
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
    k_m = 100
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
        r    = rand(1)
        fi   = 2.0*π*r[1]
        # Gaussian spectrum
        flag = true;
        while flag
            r = rand(1)
            r = r[1]        
            k = k_m*r;
            d = k*k*exp(-0.5*k*k);
            r = rand(1)
            r = r[1]
            if (r*2*exp(-1))<d
                flag = false;
            end
        end 
    
        k     = sqrt(2)*k;
        r     = rand(1)
        r     = r[1]    
        theta = acos(1-2*r);
        V1    = k*sin(fi)*sin(theta)/If[1];
        V2    = k*cos(fi)*sin(theta)/If[2];
        V3    = k*cos(theta)/If[3]; 
        r     = randn(1)
        a     = r[1]
        r     = randn(1)
        b     = r[1]
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