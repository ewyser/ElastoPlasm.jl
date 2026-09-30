
# The elastic laws live with the typed tensors (`_trial_elastic_strain`/`_trial_elastic_stress`
# in `strain.jl`/`stress.jl`). Stress and strain are full 3D tensors in every dimension; the
# kinematics (`ΔF`, `ω`) stay D×D and enter through the tensor constructors (plane strain in 2D).

"""
    elast(mpts::Point{T1,T2,D,CM,ST,EL}) where {EL<:HenckySolid}

Finite-strain elastic predictor, linear-Hencky law: push the stored logarithmic
strain forward through `ΔFᵢⱼ` and evaluate the trial Kirchhoff stress from it. Writes
a `LogarithmicStrain` into `mpts.s.ϵᵢ[p]` and a `KirchhoffStress` into
`mpts.s.τᵢ[p]`. Which of this method / the `EL<:ImprovedHenckySolid` method below
runs is picked by `Point`'s own `EL` type parameter (`solver.material.elastic`, see
`AbstractElasticLaw` in `lagrangian.jl`) — no branch anywhere, same pattern `retmap`
already uses for `Point`'s `CM` type parameter.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,ST,EL}) where {T1,T2,D,CM,ST<:LogarithmicStrain,EL<:HenckySolid}
    p = @index(Global)
    if p ≤ mpts.nmp
        cmp           = mpts.s.cmp[p]
        ϵᵢ           = _trial_elastic_strain(mpts.s.ΔFᵢⱼ[p], mpts.s.ϵᵢ[p])
        mpts.s.ϵᵢ[p] = ϵᵢ
        mpts.s.τᵢ[p] = _trial_elastic_stress(ϵᵢ, cmp)
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
        ϵᵢ           = _trial_elastic_strain(mpts.s.ΔFᵢⱼ[p], mpts.s.ϵᵢ[p])
        mpts.s.ϵᵢ[p] = ϵᵢ
        mpts.s.τᵢ[p] = _trial_elastic_stress_improved(ϵᵢ, cmp, mpts.n₀[p])
    end
end

"""
    elast(mpts::Point{T1,T2,D,CM,ST}) where {ST<:InfinitesimalStrain}

Infinitesimal (small-strain) elastic update at material points: Jaumann-rate Cauchy
stress increment `σ ← σ + Del·ϵ + (σω' + σ'ω)`, with `ϵ = ½(ΔF + ΔFᵀ) − I` and the spin `ω`,
both built on the in-plane D×D block (zero out of plane). Writes an `InfinitesimalStrain` into `mpts.s.ϵᵢ[p]` and a
`CauchyStress` into `mpts.s.σᵢ[p]`; the increment is summed on the pressure/deviator split,
where `Del·ϵ` is `-Kc·vol` on the pressure and `Gc` times the (engineering) strain deviator,
doubled on the normals.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,ST}) where {T1,T2,D,CM,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        cmp = mpts.s.cmp[p]
        ∇v  = mpts.s.∇vᵢⱼ[p]
        ω   = T2(0.5) .* (∇v - ∇v')
        # plane strain (2D) / uniaxial strain (1D): no out-of-plane stretch (ΔF_zz = 1) or spin, so
        # the strain increment and the Jaumann term are zero out of plane and are built in-plane
        ΔF  = mpts.s.ΔFᵢⱼ[p]
        ϵ   = InfinitesimalStrain(T2(0.5) .* (ΔF + ΔF') .- SMatrix{D,D,T2}(I))
        σJ  = get_tensor(mpts.s.σᵢ[p], Val(D))
        jaumann       = CauchyStress(σJ * ω' + σJ' * ω)
        mpts.s.ϵᵢ[p] = ϵ
        mpts.s.ωᵢⱼ[p] = eltype(mpts.s.ωᵢⱼ)(ω)
        σ             = mpts.s.σᵢ[p]
        mpts.s.σᵢ[p] = CauchyStress(σ.p + jaumann.p - cmp.Kc * ϵ.vol, σ.dev + jaumann.dev + cmp.Gc * ϵ.dev .* _shear_weights(T2))
    end
end
