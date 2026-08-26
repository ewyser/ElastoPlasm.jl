"""
    voigt_of(M::SMatrix{2,2,T}) -> SVector{3,T}
    voigt_of(M::SMatrix{3,3,T}) -> SVector{6,T}

Stress-convention Voigt view of an arbitrary symmetric matrix with no typed-tensor
home — e.g. a Jaumann-rate correction term, not itself "the" particle's stress (so
`get_tensor`/`get_voigt` on a stored `CauchyStress`/`KirchhoffStress` don't apply).
Replaces the old free functions `mutate`/`_mutate`, which used to also handle the
strain-side engineering↔tensor shear conversion — that conversion now lives on
`LogarithmicStrain`/`InfinitesimalStrain`'s `SVector` constructors and `get_voigt`
directly (`tensor.jl`), since `Del`-facing strain values have a real typed home this
one genuinely doesn't.
"""
@inline voigt_of(M::SMatrix{2,2,T}) where {T} = SVector{3,T}(M[1,1], M[2,2], M[1,2])
@inline voigt_of(M::SMatrix{3,3,T}) where {T} = SVector{6,T}(M[1,1], M[2,2], M[3,3], M[2,3], M[1,3], M[1,2])

# The former free functions `_logarithmic_strain`/`_kirchoff_stress` are gone: their
# math now lives on the typed tensors as `_trial_elastic_strain`/`_trial_elastic_stress`
# (`src/boot/needs/types/tensor.jl`), which is where the volumetric/deviatoric
# split naturally belongs now that it is what actually gets stored.

"""
    elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {ST<:LogarithmicStrain}

Finite-strain elastic predictor: push the stored logarithmic strain forward through
`ΔFᵢⱼ` and evaluate the trial Kirchhoff stress from it. Writes a `LogarithmicStrain`
into `mpts.s.ϵᵢⱼ[p]` and a `KirchhoffStress` into `mpts.s.τᵢⱼ[p]`.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM,TM,TV,TS,ST<:LogarithmicStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        cmp           = mpts.s.cmp[p]
        ϵᵢⱼ           = _trial_elastic_strain(mpts.s.ΔFᵢⱼ[p], mpts.s.ϵᵢⱼ[p])
        τᵢⱼ           = _trial_elastic_stress(ϵᵢⱼ, cmp.Kc, cmp.Gc)
        mpts.s.ϵᵢⱼ[p] = ϵᵢⱼ
        mpts.s.τᵢⱼ[p] = τᵢⱼ
    end
end

"""
    elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {ST<:InfinitesimalStrain}

Infinitesimal (small-strain) elastic update at material points: Jaumann-rate Cauchy
stress increment `σ ← σ + Del·ϵ + (σω' + σ'ω)`. Writes an `InfinitesimalStrain` into
`mpts.s.ϵᵢⱼ[p]` and a `CauchyStress` into `mpts.s.σᵢⱼ[p]`; the incremental arithmetic
itself still happens in Voigt `SVector` form (via `get_voigt`) so the numbers are
unchanged, with the result wrapped once at the point of the store. `Del` expects the
engineering-Voigt strain vector, which `get_voigt(InfinitesimalStrain(ϵ))` now
produces directly (see `tensor.jl`) — `ϵ` itself stays a raw tensor-shear `SMatrix`
until wrapped, same as before.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,TM,TV,TS,ST}) where {T1,T2,D,CM,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        Del = mpts.s.cmp[p].Del
        ΔF  = mpts.s.ΔFᵢⱼ[p]
        ∇v  = mpts.s.∇vᵢⱼ[p]
        ϵ   = T2(0.5) .* (ΔF + ΔF') .- eltype(mpts.s.ΔFᵢⱼ)(I)
        ω   = T2(0.5) .* (∇v - ∇v')
        σ   = get_voigt(mpts.s.σᵢⱼ[p])
        σJ  = get_tensor(mpts.s.σᵢⱼ[p])
        jaumann       = σJ * ω' + σJ' * ω
        mpts.s.ϵᵢⱼ[p] = InfinitesimalStrain(eltype(mpts.s.ΔFᵢⱼ)(ϵ))
        mpts.s.ωᵢⱼ[p] = eltype(mpts.s.ωᵢⱼ)(ω)
        mpts.s.σᵢⱼ[p]  = CauchyStress(σ + typeof(σ)(Del * get_voigt(InfinitesimalStrain(ϵ)) .+ voigt_of(jaumann)))
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
        σ = get_voigt(mpts.s.σᵢⱼ[p])
        mpts.s.σᵢⱼ[p] = CauchyStress(SVector{3,T2}(
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
        σ = get_voigt(mpts.s.σᵢⱼ[p])
        mpts.s.σᵢⱼ[p] = CauchyStress(SVector{6,T2}(
            σ[1] + (Del[1,1]*ϵxx+Del[1,2]*ϵyy+Del[1,3]*ϵzz+Del[1,4]*ϵyz+Del[1,5]*ϵxz+Del[1,6]*ϵxy) + T2(2.0)*(ωxy*σ[6]+ωxz*σ[5]),
            σ[2] + (Del[2,1]*ϵxx+Del[2,2]*ϵyy+Del[2,3]*ϵzz+Del[2,4]*ϵyz+Del[2,5]*ϵxz+Del[2,6]*ϵxy) - T2(2.0)*(ωxy*σ[6]-ωyz*σ[4]),
            σ[3] + (Del[3,1]*ϵxx+Del[3,2]*ϵyy+Del[3,3]*ϵzz+Del[3,4]*ϵyz+Del[3,5]*ϵxz+Del[3,6]*ϵxy) - T2(2.0)*(ωxz*σ[5]+ωyz*σ[4]),
            σ[4] + (Del[4,1]*ϵxx+Del[4,2]*ϵyy+Del[4,3]*ϵzz+Del[4,4]*ϵyz+Del[4,5]*ϵxz+Del[4,6]*ϵxy) + (ωyz*(σ[3]-σ[2])-σ[6]*ωxz-σ[5]*ωxy),
            σ[5] + (Del[5,1]*ϵxx+Del[5,2]*ϵyy+Del[5,3]*ϵzz+Del[5,4]*ϵyz+Del[5,5]*ϵxz+Del[5,6]*ϵxy) + (ωxz*(σ[3]-σ[1])+σ[4]*ωxy-σ[6]*ωyz),
            σ[6] + (Del[6,1]*ϵxx+Del[6,2]*ϵyy+Del[6,3]*ϵzz+Del[6,4]*ϵyz+Del[6,5]*ϵxz+Del[6,6]*ϵxy) + (ωxy*(σ[2]-σ[1])+σ[4]*ωxz+σ[5]*ωyz),
        ))
    end
end
