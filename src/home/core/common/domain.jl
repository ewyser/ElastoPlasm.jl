"""
    particle_domain(kind, F, ℓ₀, ℓ) -> ℓ

A material point's domain half-lengths once its deformation gradient `F` has changed, by basis kind: GIMP's
are its initial ones `ℓ₀` scaled by the diagonal of the stretch tensor `U = √(FᵀF)` (Charlton et al. 2017,
as AMPLE, Coombs & Augarde 2020 Eq 34); the other kinds have none to update and keep `ℓ`. Called where `F` is
updated: the explicit `deform` and the dynamic relaxation's end of a load step (`update`).
"""
@inline function particle_domain(kind::AbstractBasis, F, ℓ₀, ℓ)
    return ℓ
end
@inline function particle_domain(kind::GimpBasis, F, ℓ₀, ℓ)
    λ, n = eigen(Symmetric(F'*F))
    return diag(n*diagm(sqrt.(λ))*n') .* ℓ₀
end
