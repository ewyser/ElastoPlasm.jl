
# The Hencky push-forward and the elastic laws are in `elasticity.jl` (`_trial_elastic_strain`,
# `_trial_elastic_stress`/`_trial_elastic_stress_improved`). Stress and strain
# are full 3D tensors in every dimension; the kinematics (`ΔF`, `ω`) stay D×D and enter through the tensor constructors (plane strain in 2D).

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
`_trial_elastic_stress_improved` in `elasticity.jl`). Otherwise identical to the
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
stress increment `σ ← σ + Del·ϵ + (Δω·σ − σ·Δω)` (Borja 2013, Eq 5.46), with the strain increment
`ϵ = ½(ΔF + ΔFᵀ) − I` and the spin increment `Δω = ½(ΔF − ΔFᵀ)`, both built on the in-plane D×D
block (zero out of plane). Writes an `InfinitesimalStrain` into `mpts.s.ϵᵢ[p]` and a
`CauchyStress` into `mpts.s.σᵢ[p]`; the increment is summed on the pressure/deviator split,
where `Del·ϵ` is `-Kc·vol` on the pressure and `2Gc` times the strain deviator, its engineering
shears halved to tensor components.
"""
@kernel inbounds = true function elast(mpts::Point{T1,T2,D,CM,ST}) where {T1,T2,D,CM,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        cmp = mpts.s.cmp[p]
        ∇v  = mpts.s.∇vᵢⱼ[p]
        ω   = T2(0.5) .* (∇v - ∇v')
        # plane strain (2D) / uniaxial strain (1D): no out-of-plane stretch (ΔF_zz = 1) or spin, so
        # the strain increment and the Jaumann rotation act on the in-plane block only
        ΔF  = mpts.s.ΔFᵢⱼ[p]
        ϵ   = InfinitesimalStrain(T2(0.5) .* (ΔF + ΔF') .- SMatrix{D,D,T2}(I), zero(InfinitesimalStrain{T2}))
        # spin increment Δt·w = ½(ΔF − ΔFᵀ), since ΔF = I + Δt·∇v; Jaumann rotation Δt·(w·σ − σ·w)
        # (Borja 2013, Eq 5.46; Huang et al. 2015, Eq 32)
        Δω  = T2(0.5) .* (ΔF - ΔF')
        σ   = mpts.s.σᵢ[p]
        σJ  = get_tensor(σ, Val(D))
        σR  = CauchyStress(σJ + Δω * σJ - σJ * Δω, σ)
        mpts.s.ϵᵢ[p] = ϵ
        mpts.s.ωᵢⱼ[p] = eltype(mpts.s.ωᵢⱼ)(ω)
        mpts.s.σᵢ[p] = CauchyStress(σR.p - cmp.Kc * ϵ.vol, σR.dev + T2(2.0) * cmp.Gc * ϵ.dev .* _shear_weights(T2))
    end
end
