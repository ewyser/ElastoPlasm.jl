"""
    particle_domain(kind, F, ℓ₀, ℓ) -> ℓ

A material point's domain half-lengths once its deformation gradient `F` has changed, by basis kind: GIMP's
follow its domain update (`kind.domain`, `basis.domain`: `StretchDomain`, `UndeformedDomain`); the other kinds
have none to update and keep `ℓ`. Called where `F` is updated: the explicit `deform` and the dynamic
relaxation's end of a load step (`update`). A further domain update is one `AbstractDomain` type and one method.
"""
@inline function particle_domain(kind::AbstractBasis, F, ℓ₀, ℓ)
    return ℓ
end
@inline function particle_domain(kind::GimpBasis, F, ℓ₀, ℓ)
    return particle_domain(kind.domain, F, ℓ₀)
end
# the stretch update: the initial half-lengths scaled by the diagonal of U = √(FᵀF) (Charlton et al. 2017, as
# AMPLE, Coombs & Augarde 2020 Eq 34)
@inline function particle_domain(::StretchDomain, F, ℓ₀)
    λ, n = eigen(Symmetric(F'*F))
    return diag(n*diagm(sqrt.(λ))*n') .* ℓ₀
end
# uGIMP: the initial half-lengths (Bardenhagen & Kober 2004)
@inline function particle_domain(::UndeformedDomain, F, ℓ₀)
    return ℓ₀
end
