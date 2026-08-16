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
    if type == :tensor
        if size(ϵ) == (3,)
            return SMatrix{2,2}(ϵ[1], Χ*ϵ[3], Χ*ϵ[3], ϵ[2])
        elseif size(ϵ) == (6,)
            return SMatrix{3,3}(ϵ[1],Χ*ϵ[6],Χ*ϵ[5], Χ*ϵ[6],ϵ[2],Χ*ϵ[4], Χ*ϵ[5],Χ*ϵ[4],ϵ[3])
        end
    elseif type == :voigt
        if size(ϵ) == (2,2)
            return SVector{3}(ϵ[1,1], ϵ[2,2], Χ*ϵ[1,2])
        elseif size(ϵ) == (3,3)
            return SVector{6}(ϵ[1,1], ϵ[2,2], ϵ[3,3], Χ*ϵ[2,3], Χ*ϵ[1,3], Χ*ϵ[1,2])
        end
    end
end

@inline function _mutate(ϵ::SMatrix{2,2,T}) where {T}
    return SVector{3,T}(ϵ[1,1], ϵ[2,2], 2.0*ϵ[1,2])
end
@inline function _mutate(ϵ::SMatrix{3,3,T}) where {T}
    return SVector{6,T}(ϵ[1,1], ϵ[2,2], ϵ[3,3], 2.0*ϵ[2,3], 2.0*ϵ[1,3], 2.0*ϵ[1,2])
end

@inline function _mutate(ϵ::SVector{3,T}) where {T}
    return SMatrix{2,2,T}(ϵ[1], 0.5*ϵ[3], 0.5*ϵ[3], ϵ[2])
end
@inline function _mutate(ϵ::SVector{6,T}) where {T}
    return SMatrix{3,3,T}(ϵ[1],0.5*ϵ[6],0.5*ϵ[5], 0.5*ϵ[6],ϵ[2],0.5*ϵ[4], 0.5*ϵ[5],0.5*ϵ[4],ϵ[3])
end

@inline function _logarithmic_strain(ΔFᵢⱼ::SMatrix{D,D,T2}, ϵᵢⱼ::SMatrix{D,D,T2}) where {D,T2}
    # Compute trial left Cauchy-Green tensor
    λ,n = eigen(Symmetric(ϵᵢⱼ))
    bᵢⱼ = ΔFᵢⱼ * (n * diagm(exp.(2.0 * λ)) * n') * ΔFᵢⱼ'
    # Compute the logarithmic strain tensor using eigen decomposition
    λ, n = eigen(Symmetric(bᵢⱼ))
    # Compute the logarithmic strain tensor
    ϵᵢⱼ  = T2(0.5) * (n * diagm(log.(λ)) * n')
    # Return updated logarithmic strain tensor
    return SMatrix{D,D,T2}(ϵᵢⱼ)
end
@inline function _kirchoff_stress(ϵᵢⱼ::SMatrix{2,2,T2}, Kc::T2, Gc::T2) where {T2}
    # Calculate volumetric strain
    ϵV  = tr(ϵᵢⱼ) / T2(3.0)
    # Calculate deviatoric strain 
    ϵxx = ϵᵢⱼ[1,1] - ϵV
    ϵyy = ϵᵢⱼ[2,2] - ϵV
    ϵxy = ϵᵢⱼ[1,2]
    # Calculate pressure (positive in compression)
    P   = - T2(3.0) * Kc * ϵV
    # Calculate Krichhoff stress tensor
    τxx = T2(2.0) * Gc * ϵxx - P
    τyy = T2(2.0) * Gc * ϵyy - P
    τxy = T2(2.0) * Gc * ϵxy
    # Return Krichhoff stress in Voigt notation
    return SVector{3,T2}(τxx, τyy, τxy)
end
@inline function _kirchoff_stress(ϵᵢ::SVector{3,T2}, Kc::T2, Gc::T2) where {T2}
    # Compute volumetric strain
    ϵV  = (ϵᵢ[1]+ϵᵢ[2]) / T2(3.0)
    # Calculate pressure (positive in compression)
    P   = - T2(3.0) * Kc * ϵV
    # Calculate Krichhoff stress tensor
    τxx = T2(2.0) * Gc * (ϵᵢ[1] - ϵV) - P
    τyy = T2(2.0) * Gc * (ϵᵢ[2] - ϵV) - P
    τxy = T2(2.0) * Gc *  ϵᵢ[3]
    # Return Krichhoff stress in Voigt notation
    return SVector{3,T2}(τxx, τyy, τxy)
end
@inline function _kirchoff_stress(ϵᵢⱼ::SMatrix{3,3,T2}, Kc::T2, Gc::T2) where {T2}
    # Compute volumetric strain
    ϵV  = tr(ϵᵢⱼ) / T2(3.0)
    # Calculate pressure (positive in compression)
    P   = - T2(3.0) * Kc * ϵV
    # Calculate Krichhoff stress tensor
    τxx = T2(2.0) * Gc * (ϵᵢⱼ[1,1] - ϵV) - P
    τyy = T2(2.0) * Gc * (ϵᵢⱼ[2,2] - ϵV) - P
    τzz = T2(2.0) * Gc * (ϵᵢⱼ[3,3] - ϵV) - P
    τxy = T2(2.0) * Gc *  ϵᵢⱼ[1,2]
    τxz = T2(2.0) * Gc *  ϵᵢⱼ[1,3]
    τyz = T2(2.0) * Gc *  ϵᵢⱼ[2,3]
    # Return Krichhoff stress in Voigt notation
    return SVector{6,T2}(τxx, τyy, τzz, τxy, τxz, τyz)
end
@inline function _kirchoff_stress(ϵᵢ::SVector{6,T2}, Kc::T2, Gc::T2) where {T2}
    # Compute volumetric strain
    ϵV  = (ϵᵢ[1]+ϵᵢ[2]+ϵᵢ[3]) / T2(3.0)
    # Calculate pressure (positive in compression)
    P   = - T2(3.0) * Kc * ϵV
    # Calculate Krichhoff stress tensor
    τxx = T2(2.0) * Gc * (ϵᵢ[1] - ϵV) - P
    τyy = T2(2.0) * Gc * (ϵᵢ[2] - ϵV) - P
    τzz = T2(2.0) * Gc * (ϵᵢ[3] - ϵV) - P
    τxy = T2(2.0) * Gc *  ϵᵢ[4]
    τxz = T2(2.0) * Gc *  ϵᵢ[5]
    τyz = T2(2.0) * Gc *  ϵᵢ[6]
    # Return Krichhoff stress in Voigt notation
    return SVector{6,T2}(τxx, τyy, τzz, τxy, τxz, τyz)
end


@kernel inbounds = true function elast(mpts::Point{T1,T2,D,E,R},Kc,Gc) where {T1,T2,D,E<:FiniteElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp
        ϵᵢⱼ           = _logarithmic_strain(mpts.s.ΔFᵢⱼ[p], mpts.s.ϵᵢⱼ[p])
        τᵢ            = _kirchoff_stress(ϵᵢⱼ, Kc, Gc)
        mpts.s.ϵᵢⱼ[p] = ϵᵢⱼ
        mpts.s.τᵢ[p]  = τᵢ
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
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,E,R},Del) where {T1,T2,D,E<:LinearElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp
        ΔF = mpts.s.ΔFᵢⱼ[p]
        ∇v = mpts.s.∇vᵢⱼ[p]
        ϵ = T2(0.5) .* (ΔF + ΔF') .- eltype(mpts.s.ΔFᵢⱼ)(I)
        ω = T2(0.5) .* (∇v - ∇v')
        σJ = mutate(mpts.s.σᵢ[p], T2(1.0), :tensor)
        jaumann = σJ * ω' + σJ' * ω
        mpts.s.ϵᵢⱼ[p] = eltype(mpts.s.ϵᵢⱼ)(ϵ)
        mpts.s.ωᵢⱼ[p]  = eltype(mpts.s.ωᵢⱼ)(ω)
        mpts.s.σᵢ[p]   = mpts.s.σᵢ[p] + eltype(mpts.s.σᵢ)(Del * mutate(ϵ, T2(2.0), :voigt) .+ mutate(jaumann, T2(1.0), :voigt))
    end
end

@kernel inbounds = true function elast_fast(mpts::Point{T1,T2,2,E,R},Del) where {T1,T2,E<:LinearElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp
        ΔF = mpts.s.ΔFᵢⱼ[p]
        ω = mpts.s.ωᵢⱼ[p]
        ϵxx = ΔF[1,1] - T2(1.0)
        ϵyy = ΔF[2,2] - T2(1.0)
        ϵxy = ΔF[1,2] + ΔF[2,1]
        ωxy = T2(0.5) * (ω[1,2] - ω[2,1])
        σ = mpts.s.σᵢ[p]
        mpts.s.σᵢ[p] = SVector{3,T2}(
            σ[1] + (Del[1,1]*ϵxx+Del[1,2]*ϵyy+Del[1,3]*ϵxy) + ωxy*T2(2.0)*σ[3],
            σ[2] + (Del[2,1]*ϵxx+Del[2,2]*ϵyy+Del[2,3]*ϵxy) - ωxy*T2(2.0)*σ[3],
            σ[3] + (Del[3,1]*ϵxx+Del[3,2]*ϵyy+Del[3,3]*ϵxy) + ωxy*T2(1.0)*(σ[2]-σ[1]),
        )
    end
end
@kernel inbounds = true function elast_fast(mpts::Point{T1,T2,3,E,R},Del) where {T1,T2,E<:LinearElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp
        ΔF = mpts.s.ΔFᵢⱼ[p]
        ω = mpts.s.ωᵢⱼ[p]
        ϵxx = ΔF[1,1] - T2(1.0)
        ϵyy = ΔF[2,2] - T2(1.0)
        ϵzz = ΔF[3,3] - T2(1.0)
        ϵyz = ΔF[2,3] + ΔF[3,2]
        ϵxz = ΔF[1,3] + ΔF[3,1]
        ϵxy = ΔF[1,2] + ΔF[2,1]
        ωyz = T2(0.5) * (ω[2,3] - ω[3,2])
        ωxz = T2(0.5) * (ω[1,3] - ω[3,1])
        ωxy = T2(0.5) * (ω[1,2] - ω[2,1])
        σ = mpts.s.σᵢ[p]
        mpts.s.σᵢ[p] = SVector{6,T2}(
            σ[1] + (Del[1,1]*ϵxx+Del[1,2]*ϵyy+Del[1,3]*ϵzz+Del[1,4]*ϵyz+Del[1,5]*ϵxz+Del[1,6]*ϵxy) + T2(2.0)*(ωxy*σ[6]+ωxz*σ[5]),
            σ[2] + (Del[2,1]*ϵxx+Del[2,2]*ϵyy+Del[2,3]*ϵzz+Del[2,4]*ϵyz+Del[2,5]*ϵxz+Del[2,6]*ϵxy) - T2(2.0)*(ωxy*σ[6]-ωyz*σ[4]),
            σ[3] + (Del[3,1]*ϵxx+Del[3,2]*ϵyy+Del[3,3]*ϵzz+Del[3,4]*ϵyz+Del[3,5]*ϵxz+Del[3,6]*ϵxy) - T2(2.0)*(ωxz*σ[5]+ωyz*σ[4]),
            σ[4] + (Del[4,1]*ϵxx+Del[4,2]*ϵyy+Del[4,3]*ϵzz+Del[4,4]*ϵyz+Del[4,5]*ϵxz+Del[4,6]*ϵxy) + (ωyz*(σ[3]-σ[2])-σ[6]*ωxz-σ[5]*ωxy),
            σ[5] + (Del[5,1]*ϵxx+Del[5,2]*ϵyy+Del[5,3]*ϵzz+Del[5,4]*ϵyz+Del[5,5]*ϵxz+Del[5,6]*ϵxy) + (ωxz*(σ[3]-σ[1])+σ[4]*ωxy-σ[6]*ωyz),
            σ[6] + (Del[6,1]*ϵxx+Del[6,2]*ϵyy+Del[6,3]*ϵzz+Del[6,4]*ϵyz+Del[6,5]*ϵxz+Del[6,6]*ϵxy) + (ωxy*(σ[2]-σ[1])+σ[4]*ωxz+σ[5]*ωyz),
        )
    end
end