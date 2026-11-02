abstract type AbstractConstitutiveModel{T2,D} end

export DruckerPrager,VonMises

"""
    DruckerPrager{T2,D} <: AbstractConstitutiveModel{T2,D}

A point's constants for the Drucker–Prager model (`material.plastic = "DP"`): its elastic moduli and yield surface,
fixed once built. The evolving plastic state (`Δλ`, `ϵpII`, `ϵpV`) is the solid phase's (`PointSolidPhase`).

# Fields
- `Gc::T2`, `Kc::T2`: the shear and bulk moduli [Pa]; the elastic laws apply them to the stress's pressure and
  deviator directly (`get_Del` builds the 6×6 Voigt stiffness when a matrix is needed).
- `Hp::T2`: the softening modulus [Pa]: the cohesion is `max(cᵣ, c₀ + Hp·ϵpII)`.
- `c₀::T2`, `cᵣ::T2`: the initial and residual cohesion [Pa].
- `ϕ₀::T2`: the friction angle [rad].

`T2` is the float type, `D` the dimension: in 2D the cone matches Mohr–Coulomb in plane strain, in 3D the inner cone.

The return map is `retmap`'s method for `Point`'s `CM <: DruckerPrager` and its strain type `ST`; the elastic law is
`Point`'s `EL`, not part of these constants.

# Constructors
    DruckerPrager{T2,D}(Gc, Kc, Hp, c₀, cᵣ, ϕ₀)

`setup_cmp` builds one per point.

See also: `VonMises`, `PLASTIC_LAWS`.
"""
struct DruckerPrager{T2,D} <: AbstractConstitutiveModel{T2,D}
    Gc ::T2
    Kc ::T2
    Hp ::T2
    c₀ ::T2
    cᵣ ::T2
    ϕ₀ ::T2
end

"""
    VonMises{T2,D} <: AbstractConstitutiveModel{T2,D}

A point's constants for the von Mises (J2) model (`material.plastic = "VM"`): `DruckerPrager`'s without the friction
angle, since the yield surface does not depend on the pressure.

# Fields
- `Gc::T2`, `Kc::T2`: the shear and bulk moduli [Pa].
- `Hp::T2`: the softening modulus [Pa]: the yield radius is `κ = max(cᵣ, c₀ + Hp·ϵpII)`.
- `c₀::T2`, `cᵣ::T2`: the initial and residual yield radius [Pa], read as a Lode radius, `√(2J₂) = κ`.

`T2` is the float type, `D` the dimension.

# Constructors
    VonMises{T2,D}(Gc, Kc, Hp, c₀, cᵣ)

`setup_cmp` builds one per point.

See also: `DruckerPrager`, `PLASTIC_LAWS`.
"""
struct VonMises{T2,D} <: AbstractConstitutiveModel{T2,D}
    Gc ::T2
    Kc ::T2
    Hp ::T2
    c₀ ::T2
    cᵣ ::T2
end
