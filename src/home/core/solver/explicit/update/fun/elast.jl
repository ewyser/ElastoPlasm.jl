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
    return SVector{3,T}(ϵ[1,1], ϵ[2,2], T(2.0)*ϵ[1,2])
end
@inline function _mutate(ϵ::SMatrix{3,3,T}) where {T}
    return SVector{6,T}(ϵ[1,1], ϵ[2,2], ϵ[3,3], T(2.0)*ϵ[2,3], T(2.0)*ϵ[1,3], T(2.0)*ϵ[1,2])
end

@inline function _mutate(ϵ::SVector{3,T}) where {T}
    return SMatrix{2,2,T}(ϵ[1], T(0.5)*ϵ[3], T(0.5)*ϵ[3], ϵ[2])
end
@inline function _mutate(ϵ::SVector{6,T}) where {T}
    return SMatrix{3,3,T}(ϵ[1],T(0.5)*ϵ[6],T(0.5)*ϵ[5], T(0.5)*ϵ[6],ϵ[2],T(0.5)*ϵ[4], T(0.5)*ϵ[5],T(0.5)*ϵ[4],ϵ[3])
end

# The former free functions `_logarithmic_strain`/`_kirchoff_stress` are gone: their
# math now lives on the typed tensors as `_trial_elastic_strain`/`_trial_elastic_stress`
# (`src/boot/needs/types/concrete/tensor.jl`), which is where the volumetric/deviatoric
# split naturally belongs now that it is what actually gets stored.

"""
    elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {ST<:LogarithmicStrain}

Finite-strain elastic predictor: push the stored logarithmic strain forward through
`ΔFᵢⱼ` and evaluate the trial Kirchhoff stress from it. Writes a `LogarithmicStrain`
into `mpts.s.ϵᵢⱼ[p]` and a `KirchhoffStress` into `mpts.s.τᵢ[p]`.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM,TM,TV,TS,ST<:LogarithmicStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        cmp           = mpts.s.cmp[p]
        ϵᵢⱼ           = _trial_elastic_strain(mpts.s.ΔFᵢⱼ[p], mpts.s.ϵᵢⱼ[p])
        τᵢ            = _trial_elastic_stress(ϵᵢⱼ, cmp.Kc, cmp.Gc)
        mpts.s.ϵᵢⱼ[p] = ϵᵢⱼ
        mpts.s.τᵢ[p]  = τᵢ
    end
end

"""
    elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {ST<:InfinitesimalStrain}

Infinitesimal (small-strain) elastic update at material points: Jaumann-rate Cauchy
stress increment `σ ← σ + Del·ϵ + (σω' + σ'ω)`. Writes an `InfinitesimalStrain` into
`mpts.s.ϵᵢⱼ[p]` and a `CauchyStress` into `mpts.s.σᵢ[p]`; the incremental arithmetic
itself still happens in Voigt `SVector` form (via `get_vector`) so the numbers are
unchanged, with the result wrapped once at the point of the store.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        Del = mpts.s.cmp[p].Del
        ΔF  = mpts.s.ΔFᵢⱼ[p]
        ∇v  = mpts.s.∇vᵢⱼ[p]
        ϵ   = T2(0.5) .* (ΔF + ΔF') .- eltype(mpts.s.ΔFᵢⱼ)(I)
        ω   = T2(0.5) .* (∇v - ∇v')
        σ   = get_vector(mpts.s.σᵢ[p])
        σJ  = mutate(σ, T2(1.0), :tensor)
        jaumann       = σJ * ω' + σJ' * ω
        mpts.s.ϵᵢⱼ[p] = InfinitesimalStrain(eltype(mpts.s.ΔFᵢⱼ)(ϵ))
        mpts.s.ωᵢⱼ[p] = eltype(mpts.s.ωᵢⱼ)(ω)
        mpts.s.σᵢ[p]  = CauchyStress(σ + typeof(σ)(Del * mutate(ϵ, T2(2.0), :voigt) .+ mutate(jaumann, T2(1.0), :voigt)))
    end
end

@kernel inbounds = true function elast_fast(mpts::Point{T1,T2,2,CM,TM,TV,TS,ST}) where {T1,T2,CM,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        Del = mpts.s.cmp[p].Del
        ΔF = mpts.s.ΔFᵢⱼ[p]
        ω = mpts.s.ωᵢⱼ[p]
        ϵxx = ΔF[1,1] - T2(1.0)
        ϵyy = ΔF[2,2] - T2(1.0)
        ϵxy = ΔF[1,2] + ΔF[2,1]
        ωxy = T2(0.5) * (ω[1,2] - ω[2,1])
        σ = get_vector(mpts.s.σᵢ[p])
        mpts.s.σᵢ[p] = CauchyStress(SVector{3,T2}(
            σ[1] + (Del[1,1]*ϵxx+Del[1,2]*ϵyy+Del[1,3]*ϵxy) + ωxy*T2(2.0)*σ[3],
            σ[2] + (Del[2,1]*ϵxx+Del[2,2]*ϵyy+Del[2,3]*ϵxy) - ωxy*T2(2.0)*σ[3],
            σ[3] + (Del[3,1]*ϵxx+Del[3,2]*ϵyy+Del[3,3]*ϵxy) + ωxy*T2(1.0)*(σ[2]-σ[1]),
        ))
    end
end
@kernel inbounds = true function elast_fast(mpts::Point{T1,T2,3,CM,TM,TV,TS,ST}) where {T1,T2,CM,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        Del = mpts.s.cmp[p].Del
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
        σ = get_vector(mpts.s.σᵢ[p])
        mpts.s.σᵢ[p] = CauchyStress(SVector{6,T2}(
            σ[1] + (Del[1,1]*ϵxx+Del[1,2]*ϵyy+Del[1,3]*ϵzz+Del[1,4]*ϵyz+Del[1,5]*ϵxz+Del[1,6]*ϵxy) + T2(2.0)*(ωxy*σ[6]+ωxz*σ[5]),
            σ[2] + (Del[2,1]*ϵxx+Del[2,2]*ϵyy+Del[2,3]*ϵzz+Del[2,4]*ϵyz+Del[2,5]*ϵxz+Del[2,6]*ϵxy) - T2(2.0)*(ωxy*σ[6]-ωyz*σ[4]),
            σ[3] + (Del[3,1]*ϵxx+Del[3,2]*ϵyy+Del[3,3]*ϵzz+Del[3,4]*ϵyz+Del[3,5]*ϵxz+Del[3,6]*ϵxy) - T2(2.0)*(ωxz*σ[5]+ωyz*σ[4]),
            σ[4] + (Del[4,1]*ϵxx+Del[4,2]*ϵyy+Del[4,3]*ϵzz+Del[4,4]*ϵyz+Del[4,5]*ϵxz+Del[4,6]*ϵxy) + (ωyz*(σ[3]-σ[2])-σ[6]*ωxz-σ[5]*ωxy),
            σ[5] + (Del[5,1]*ϵxx+Del[5,2]*ϵyy+Del[5,3]*ϵzz+Del[5,4]*ϵyz+Del[5,5]*ϵxz+Del[5,6]*ϵxy) + (ωxz*(σ[3]-σ[1])+σ[4]*ωxy-σ[6]*ωyz),
            σ[6] + (Del[6,1]*ϵxx+Del[6,2]*ϵyy+Del[6,3]*ϵzz+Del[6,4]*ϵyz+Del[6,5]*ϵxz+Del[6,6]*ϵxy) + (ωxy*(σ[2]-σ[1])+σ[4]*ωxz+σ[5]*ωyz),
        ))
    end
end
