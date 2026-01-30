"""
    mutate(ϵ, Χ, type)

Convert between tensor and Voigt notation for strain or stress, or mutate tensor forms.

# Arguments
- `ϵ`: Strain or stress tensor or vector.
- `Χ`: Scaling factor (e.g., 1.0 or 2.0).
- `type`: Either `:tensor` or `:voigt` for conversion type.

# Returns
- `ϵmut`: Mutated tensor or vector in the requested form.
"""
@views function mutate(ϵ,Χ,type)
    if type == :tensor # α = 1/2 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (3,)
            ϵmut = [  ϵ[1] Χ*ϵ[3];
                    Χ*ϵ[3]   ϵ[2]]
        elseif size(ϵ) == (6,)
            ϵmut = [  ϵ[1] Χ*ϵ[6] Χ*ϵ[5];
                    Χ*ϵ[6]   ϵ[2] Χ*ϵ[4];
                    Χ*ϵ[5] Χ*ϵ[4]   ϵ[3]]
        end
    elseif type == :voigt # α = 2.0 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (2,2)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],Χ*ϵ[1,2]) #xx,yy,zz,xy
        elseif size(ϵ) == (3,3)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],ϵ[3,3],Χ*ϵ[2,3],Χ*ϵ[1,3],Χ*ϵ[1,2]) #xx,yy,zz,yz,xz,xy
        end
    end
    return ϵmut
end
"""
    finite_elast(mpts::Point{T1,T2}, Del) where {T1,T2}

Kernel for finite deformation elasticity update at material points.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `Del`: Elasticity matrix.

# Returns
- Updates stress and strain fields in-place.
"""
@views @kernel inbounds = true function elast(mpts::Point{T1,T2,D,E,R},Del) where {T1,T2,D,E<:FiniteElasticity,R}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # update left cauchy-green tensor
        mpts.s.bᵢⱼ[:,:,p].= mpts.s.ΔFᵢⱼ[:,:,p]*mpts.s.bᵢⱼ[:,:,p]*mpts.s.ΔFᵢⱼ[:,:,p]'
        # compute logarithmic strain tensor
        λ,n               = eigen(mpts.s.bᵢⱼ[:,:,p],sortby=nothing)
        mpts.s.ϵᵢⱼ[:,:,p].= T2(0.5).*(n*diagm(log.(λ))*n')
        # krichhoff stress tensor
        mpts.s.τᵢ[:,p]    = Del*mutate(mpts.s.ϵᵢⱼ[:,:,p],T2(2.0),:voigt)
    end
end
"""
    infinitesimal_elast(mpts::Point{T1,T2}, Del) where {T1,T2}

Kernel for infinitesimal (small strain) elasticity update at material points.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `Del`: Elasticity matrix.

# Returns
- Updates stress and strain fields in-place.
"""
@views @kernel inbounds = true function elast(mpts::Point{T1,T2,D,E,R},Del) where {T1,T2,D,E<:LinearElasticity,R}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # calculate elastic strains & spins
        mpts.s.ϵᵢⱼ[:,:,p] .= T2(0.5).*(mpts.s.ΔFᵢⱼ[:,:,p]+mpts.s.ΔFᵢⱼ[:,:,p]').-mpts.δᵢⱼ[:,:] 
        mpts.s.ωᵢⱼ[:,:,p] .= T2(0.5).*(mpts.s.∇vᵢⱼ[:,:,p]-mpts.s.∇vᵢⱼ[:,:,p]')
        # update cauchy stress tensor
        mpts.s.σJᵢⱼ[:,:,p].= mutate(mpts.s.σᵢ[:,p],T2(1.0),:tensor)
        mpts.s.σJᵢⱼ[:,:,p].= mpts.s.σJᵢⱼ[:,:,p]*mpts.s.ωᵢⱼ[:,:,p]'+mpts.s.σJᵢⱼ[:,:,p]'*mpts.s.ωᵢⱼ[:,:,p]
        mpts.s.σᵢ[:,p]   .+= Del*mutate(mpts.s.ϵᵢⱼ[:,:,p],T2(2.0),:voigt).+mutate(mpts.s.σJᵢⱼ[:,:,p],T2(1.0),:voigt)
    end  
end

@kernel inbounds = true function elast_fast(mpts::Point{T1,T2,D,E,R},Del) where {T1,T2,D<:TwoDimension,E<:LinearElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # calculate elastic strains
        ϵxx  = (mpts.s.ΔFᵢⱼ[1,1,p])-T2(1.0)
        ϵyy  = (mpts.s.ΔFᵢⱼ[2,2,p])-T2(1.0)
        ϵxy  = (mpts.s.ΔFᵢⱼ[1,2,p]+mpts.s.ΔFᵢⱼ[2,1,p])
        # calculate elastic spins
        ωxy  = T2(0.5)*(mpts.s.ωᵢⱼ[1,2,p]-mpts.s.ωᵢⱼ[2,1,p]) 
        # update jaumann stress tensor
        σxx0 = copy(mpts.s.σᵢ[1,p])
        σyy0 = copy(mpts.s.σᵢ[2,p])
        σxy0 = copy(mpts.s.σᵢ[3,p])
        # update cauchy stress tensor
        mpts.s.σᵢ[1,p]+= (Del[1,1]*ϵxx+Del[1,2]*ϵyy+Del[1,3]*ϵxy)+ωxy*T2(2.0)*σxy0
        mpts.s.σᵢ[2,p]+= (Del[2,1]*ϵxx+Del[2,2]*ϵyy+Del[2,3]*ϵxy)-ωxy*T2(2.0)*σxy0
        mpts.s.σᵢ[3,p]+= (Del[3,1]*ϵxx+Del[3,2]*ϵyy+Del[3,3]*ϵxy)+ωxy*T2(1.0)*(σyy0-σxx0)
    end  
end
@kernel inbounds = true function elast_fast(mpts::Point{T1,T2,D,E,R},Del) where {T1,T2,D<:ThreeDimension,E<:LinearElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # calculate elastic strains
        ϵxx  = (mpts.s.ΔFᵢⱼ[1,1,p])-T2(1.0)
        ϵyy  = (mpts.s.ΔFᵢⱼ[2,2,p])-T2(1.0)
        ϵzz  = (mpts.s.ΔFᵢⱼ[3,3,p])-T2(1.0)      
        ϵyz  = (mpts.s.ΔFᵢⱼ[2,3,p]+mpts.s.ΔFᵢⱼ[3,2,p])
        ϵxz  = (mpts.s.ΔFᵢⱼ[1,3,p]+mpts.s.ΔFᵢⱼ[3,1,p])
        ϵxy  = (mpts.s.ΔFᵢⱼ[1,2,p]+mpts.s.ΔFᵢⱼ[2,1,p])
        # calculate elastic spins
        ωyz  = T2(0.5)*(mpts.s.ωᵢⱼ[2,3,p]-mpts.s.ωᵢⱼ[3,2,p]) 
        ωxz  = T2(0.5)*(mpts.s.ωᵢⱼ[1,3,p]-mpts.s.ωᵢⱼ[3,1,p]) 
        ωxy  = T2(0.5)*(mpts.s.ωᵢⱼ[1,2,p]-mpts.s.ωᵢⱼ[2,1,p]) 
        # update jaumann stress tensor
        σxx0 = copy(mpts.s.σᵢ[1,p])
        σyy0 = copy(mpts.s.σᵢ[2,p])
        σzz0 = copy(mpts.s.σᵢ[3,p])
        σyz0 = copy(mpts.s.σᵢ[4,p])
        σxz0 = copy(mpts.s.σᵢ[5,p])
        σxy0 = copy(mpts.s.σᵢ[6,p])
        
        #=
        mpts.s.σᵢ[1,p]+= T2(2.0)*(ωxy*σxy0 + ωxz*σxz0)
        mpts.s.σᵢ[2,p]-= T2(2.0)*(ωxy*σxy0 - ωyz*σyz0)
        mpts.s.σᵢ[3,p]-= T2(2.0)*(ωxz*σxz0 + ωyz*σyz0)
        
        mpts.s.σᵢ[6,p]+= ( ωxy*(σyy0-σxx0)+σyz0*ωxz+σxz0*ωyz ) # xy
        mpts.s.σᵢ[4,p]+= ( ωyz*(σzz0-σyy0)-σxy0*ωxz-σxz0*ωxy ) # yz 
        mpts.s.σᵢ[5,p]+= ( ωxz*(σzz0-σxx0)+σyz0*ωxy-σxy0*ωyz ) # xz
        # update cauchy stress tensor
        # σxx, σyy, σzz = mpts.s.σᵢ[1,p]       , mpts.s.σᵢ[2,p]  , mpts.s.σᵢ[3,p]
        # σyx, σzy, σzx = mpts.s.σᵢ[6,p]       , mpts.s.σᵢ[4,p]  , mpts.s.σᵢ[5,p]
# 4 zy
# 5 zx
# 6 xy

        mpts.s.σᵢ[1,p]+= (Del[1,1]*ϵxx+Del[1,2]*ϵyy+Del[1,3]*ϵzz+Del[1,4]*ϵyz+Del[1,5]*ϵxz+Del[1,6]*ϵxy)
        mpts.s.σᵢ[2,p]+= (Del[2,1]*ϵxx+Del[2,2]*ϵyy+Del[2,3]*ϵzz+Del[2,4]*ϵyz+Del[1,5]*ϵxz+Del[1,6]*ϵxy)
        mpts.s.σᵢ[3,p]+= (Del[3,1]*ϵxx+Del[3,2]*ϵyy+Del[3,3]*ϵzz+Del[3,4]*ϵyz+Del[1,5]*ϵxz+Del[1,6]*ϵxy)
=#
    end  
end