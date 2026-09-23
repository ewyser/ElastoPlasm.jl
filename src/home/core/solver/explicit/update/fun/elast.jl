"""
    voigt_of(M::SMatrix{2,2,T}) -> SVector{3,T}
    voigt_of(M::SMatrix{3,3,T}) -> SVector{6,T}

Stress-convention Voigt view of an arbitrary symmetric matrix with no typed-tensor
home — e.g. a Jaumann-rate correction term, not itself "the" particle's stress (so
`get_tensor`/`get_voigt` on a stored `CauchyStress`/`KirchhoffStress` don't apply).
Replaces the old free functions `mutate`/`_mutate`, which used to also handle the
strain-side engineering↔tensor shear conversion — that conversion now lives on
`LogarithmicStrain`/`InfinitesimalStrain`'s `SVector` constructors and `get_voigt`
directly (`strain.jl`), since `Del`-facing strain values have a real typed home this
one genuinely doesn't.
"""
@inline voigt_of(M::SMatrix{2,2,T}) where {T} = SVector{3,T}(M[1,1], M[2,2], M[1,2])
@inline voigt_of(M::SMatrix{3,3,T}) where {T} = SVector{6,T}(M[1,1], M[2,2], M[3,3], M[2,3], M[1,3], M[1,2])

# The former free functions `_logarithmic_strain`/`_kirchoff_stress` are gone: their
# math now lives on the typed tensors as `_trial_elastic_strain`/`_trial_elastic_stress`
# (`src/boot/needs/types/problem/strain.jl`/`stress.jl`), which is where the
# volumetric/deviatoric split naturally belongs now that it is what actually gets stored.

"""
    elast(mpts::Point{T1,T2,D,CM,ST,EL}) where {EL<:HenckySolid}

Finite-strain elastic predictor, linear-Hencky law: push the stored logarithmic
strain forward through `ΔFᵢⱼ` and evaluate the trial Kirchhoff stress from it. Writes
a `LogarithmicStrain` into `mpts.s.ϵᵢⱼ[p]` and a `KirchhoffStress` into
`mpts.s.τᵢⱼ[p]`. Which of this method / the `EL<:ImprovedHenckySolid` method below
runs is picked by `Point`'s own `EL` type parameter (`solver.material.elastic`, see
`AbstractElasticLaw` in `lagrangian.jl`) — no branch anywhere, same pattern `retmap`
already uses for `Point`'s `CM` type parameter.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,ST,EL}) where {T1,T2,D,CM,ST<:LogarithmicStrain,EL<:HenckySolid}
    p = @index(Global)
    if p ≤ mpts.nmp
        cmp           = mpts.s.cmp[p]
        ϵᵢⱼ           = _trial_elastic_strain(mpts.s.ΔFᵢⱼ[p], mpts.s.ϵᵢⱼ[p])
        τᵢⱼ           = _trial_elastic_stress(ϵᵢⱼ, cmp)
        mpts.s.ϵᵢⱼ[p] = ϵᵢⱼ
        mpts.s.τᵢⱼ[p] = τᵢⱼ
    end
end

"""
    elast(mpts::Point{T1,T2,D,CM,ST,EL}) where {EL<:ImprovedHenckySolid}

Finite-strain elastic predictor, porosity-weighted "improved Hencky" law (see
`_trial_elastic_stress_improved` in `stress.jl`). Otherwise identical to the
`EL<:HenckySolid` method above.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,ST,EL}) where {T1,T2,D,CM,ST<:LogarithmicStrain,EL<:ImprovedHenckySolid}
    p = @index(Global)
    if p ≤ mpts.nmp
        cmp           = mpts.s.cmp[p]
        ϵᵢⱼ           = _trial_elastic_strain(mpts.s.ΔFᵢⱼ[p], mpts.s.ϵᵢⱼ[p])
        τᵢⱼ           = _trial_elastic_stress_improved(ϵᵢⱼ, cmp, mpts.n₀[p])
        mpts.s.ϵᵢⱼ[p] = ϵᵢⱼ
        mpts.s.τᵢⱼ[p] = τᵢⱼ
    end
end

"""
    elast(mpts::Point{T1,T2,D,CM,ST}) where {ST<:InfinitesimalStrain}

Infinitesimal (small-strain) elastic update at material points: Jaumann-rate Cauchy
stress increment `σ ← σ + Del·ϵ + (σω' + σ'ω)`. Writes an `InfinitesimalStrain` into
`mpts.s.ϵᵢⱼ[p]` and a `CauchyStress` into `mpts.s.σᵢⱼ[p]`; the incremental arithmetic
itself still happens in Voigt `SVector` form (via `get_voigt`) so the numbers are
unchanged, with the result wrapped once at the point of the store. `Del` expects the
engineering-Voigt strain vector, which `get_voigt(InfinitesimalStrain(ϵ))` now
produces directly (see `strain.jl`) — `ϵ` itself stays a raw tensor-shear `SMatrix`
until wrapped, same as before.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,ST}) where {T1,T2,D,CM,ST<:InfinitesimalStrain}
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
