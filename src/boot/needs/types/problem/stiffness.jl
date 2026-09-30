# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Isotropic elastic stiffness in Voigt space
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export StressStrainStiffness, get_Del

"""
    StressStrainStiffness{T}()

Volumetric and deviatoric projections on a 6-component Voigt strain `(xx, yy, zz, γyz, γxz, γxy)`
with engineering shear (`γ = 2ϵ`): `vol·ϵ` is `tr(ϵ)` on each normal component, and `dev·ϵ` is the
deviatoric strain, with `½` on the shear diagonal so that `2G·dev` maps `γ` to the shear stress
`G·γ`. Any isotropic stiffness is `K·vol + 2G·dev`; see `get_Del`.
"""
struct StressStrainStiffness{T}
    vol::SMatrix{6,6,T,36}
    dev::SMatrix{6,6,T,36}
end
@inline function StressStrainStiffness{T}() where {T}
    a, b, h = T(2.0/3.0), T(-1.0/3.0), T(0.5)
    vol = SMatrix{6,6,T,36}(
        1, 1, 1, 0, 0, 0,
        1, 1, 1, 0, 0, 0,
        1, 1, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
    )
    dev = SMatrix{6,6,T,36}(
        a, b, b, 0, 0, 0,
        b, a, b, 0, 0, 0,
        b, b, a, 0, 0, 0,
        0, 0, 0, h, 0, 0,
        0, 0, 0, 0, h, 0,
        0, 0, 0, 0, 0, h,
    )
    return StressStrainStiffness{T}(vol, dev)
end

"""
    get_Del(cmp::AbstractConstitutiveModel) -> SMatrix{6,6}

Isotropic elastic stiffness `Del = Kc·vol + 2Gc·dev` of a constitutive model, built on the fly from
its moduli (acts on engineering-shear Voigt strain, gives Voigt stress). Not stored per particle.
"""
@inline function get_Del(cmp::AbstractConstitutiveModel{T}) where {T}
    P = StressStrainStiffness{T}()
    return cmp.Kc * P.vol + T(2.0) * cmp.Gc * P.dev
end
