# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Strain/stress tensor types and subtypes
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractTensor, AbstractStrain, AbstractStress
export InfinitesimalStrain, LogarithmicStrain, CauchyStress, KirchhoffStress
export get_tensor, get_voigt, get_J2, get_τII

"""
    AbstractTensor{T}

Supertype of every per-particle second-order tensor stored on `PointSolidPhase`
(`ϵᵢⱼ`/`ϵn`/`σᵢ`/`σn`/`τᵢ`). Stored as a volumetric+deviatoric split — `(vol, dev)` for
strains, `(p, dev)` for stresses — so the physical meaning is carried by the type.

`dev` is **not** required to be trace-free: `_trial_elastic_stress` uses `tr(ϵ)/3` even
in 2D (plane-strain), while the Voigt constructors below use `-tr(σ)/S`. Always read a
tensor via `get_voigt`/`get_tensor`, never `.dev`/`.p` directly — `get_J2`/`get_τII`
re-derive a trace-free deviator from `get_voigt` for exactly this reason.

`get_voigt` is one function dispatched onto `AbstractStrain`/`AbstractStress`: strain
uses the engineering-shear convention (`γxy = 2ϵxy`), stress does not.
"""
AbstractTensor

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Strain tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    AbstractStrain{S,T,L} <: AbstractTensor{T}

Second-order strain tensor stored as `vol::T` (volumetric part) + `dev::SMatrix{S,S,T,L}`
(deviatoric part). `S` is the spatial dimension, `L == S*S`. Declared in
`types/abstract.jl`; concrete subtypes and methods live here.
"""
AbstractStrain

Base.getindex(strain::AbstractStrain{S,T,L}, i::Int, j::Int) where {S,T,L} =
    strain.dev[i,j] + (i == j ? strain.vol : zero(T))

"""
    get_tensor(strain::AbstractStrain{S,T,L}) -> SMatrix{S,S,T,L}

Reassemble the full strain tensor `dev + vol·I`.
"""
@inline function get_tensor(strain::AbstractStrain{S,T,L}) where {S,T,L}
    return strain.dev + strain.vol * SMatrix{S,S,T,L}(I)
end

"""
    get_voigt(strain::AbstractStrain) -> SVector{3,T} (2D) / SVector{6,T} (3D)

Voigt strain vector `(ϵxx,ϵyy,γxy)` / `(ϵxx,ϵyy,ϵzz,γyz,γxz,γxy)`, engineering-shear
convention (`γxy = 2ϵxy`) — required for `Del`-facing algebra (`σ = Del·ε_voigt`).
Round-trips with `LogarithmicStrain(ϵ::SVector)`/`InfinitesimalStrain(ϵ::SVector)`.
"""
@inline function get_voigt(strain::AbstractStrain{2,T,L}) where {T,L}
    return SVector{3,T}(
        strain.dev[1,1] + strain.vol,
        strain.dev[2,2] + strain.vol,
        T(2.0) * strain.dev[1,2],
    )
end
@inline function get_voigt(strain::AbstractStrain{3,T,L}) where {T,L}
    return SVector{6,T}(
        strain.dev[1,1] + strain.vol,
        strain.dev[2,2] + strain.vol,
        strain.dev[3,3] + strain.vol,
        T(2.0) * strain.dev[2,3],
        T(2.0) * strain.dev[1,3],
        T(2.0) * strain.dev[1,2],
    )
end

"""
    InfinitesimalStrain{S,T,L} <: AbstractStrain{S,T,L}

Small-strain tensor `ϵ = ½(ΔF + ΔFᵀ) - I`, used under `strain.deform == "infinitesimal"`.
Built by `_infinitesimal_strain`.
"""
struct InfinitesimalStrain{S,T,L} <: AbstractStrain{S,T,L}
    vol::T
    dev::SMatrix{S,S,T,L}
    InfinitesimalStrain{S,T,L}(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(vol, dev)
    InfinitesimalStrain(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L} = new{S,T,L}(vol, dev)
end
@adapt_struct InfinitesimalStrain

"""
    LogarithmicStrain{S,T,L} <: AbstractStrain{S,T,L}

Logarithmic (Hencky) elastic strain tensor, used under `strain.deform == "finite"`.
Built by `_trial_elastic_strain`.
"""
struct LogarithmicStrain{S,T,L} <: AbstractStrain{S,T,L}
    vol::T
    dev::SMatrix{S,S,T,L}
    function LogarithmicStrain{S,T,L}(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(vol, dev)
    end
    function LogarithmicStrain(vol::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(vol, dev)
    end
end
@adapt_struct LogarithmicStrain

# dev given as a Voigt deviator only (no volumetric part)
@inline function LogarithmicStrain(vol::T, dev::SVector{3,T}) where {T}
    return LogarithmicStrain(vol, SMatrix{2,2,T,4}(
        dev[1], dev[3],
        dev[3], dev[2],
    ))
end
# 3D: SMatrix fills column-major — dev = [dev1 dev4 dev5; dev4 dev2 dev6; dev5 dev6 dev3].
@inline function LogarithmicStrain(vol::T, dev::SVector{6,T}) where {T}
    return LogarithmicStrain(vol, SMatrix{3,3,T,9}(
        dev[1], dev[4], dev[5],
        dev[4], dev[2], dev[6],
        dev[5], dev[6], dev[3],
    ))
end

"""
    LogarithmicStrain(ϵᵢⱼ::SMatrix{S,S,T,L})

Split a full logarithmic strain tensor into `vol = tr(ϵ)/3` + `dev = ϵ - vol·I`.
"""
@inline function LogarithmicStrain(ϵᵢⱼ::SMatrix{S,S,T,L}) where {S,T,L}
    vol = tr(ϵᵢⱼ) / T(3.0)
    return LogarithmicStrain(vol, SMatrix{S,S,T,L}(ϵᵢⱼ - vol * SMatrix{S,S,T,L}(I)))
end

"""
    LogarithmicStrain(ϵ::SVector{3,T}) / LogarithmicStrain(ϵ::SVector{6,T})

Build directly from a full engineering-Voigt strain vector. `vol = tr(ϵ)/3` (same
plane-strain convention as the `SMatrix` constructor above — not stress's `tr/S`).
Exact inverse of `get_voigt`.
"""
@inline function LogarithmicStrain(ϵ::SVector{3,T}) where {T}
    vol = (ϵ[1] + ϵ[2]) / T(3.0)
    return LogarithmicStrain(vol, SMatrix{2,2,T,4}(ϵ[1] - vol, ϵ[3]/T(2.0), ϵ[3]/T(2.0), ϵ[2] - vol))
end
@inline function LogarithmicStrain(ϵ::SVector{6,T}) where {T}
    vol = (ϵ[1] + ϵ[2] + ϵ[3]) / T(3.0)
    return LogarithmicStrain(vol, SMatrix{3,3,T,9}(
        ϵ[1] - vol, ϵ[6]/T(2.0), ϵ[5]/T(2.0),
        ϵ[6]/T(2.0), ϵ[2] - vol, ϵ[4]/T(2.0),
        ϵ[5]/T(2.0), ϵ[4]/T(2.0), ϵ[3] - vol,
    ))
end

Base.zero(::Type{InfinitesimalStrain{S,T,L}}) where {S,T,L} = InfinitesimalStrain(zero(T), zero(SMatrix{S,S,T,L}))
Base.zero(::Type{LogarithmicStrain{S,T,L}})   where {S,T,L} = LogarithmicStrain(  zero(T), zero(SMatrix{S,S,T,L}))

"""
    LinearAlgebra.eigen(strain::AbstractStrain)

Eigen-decomposition of the reassembled (symmetric) strain tensor. Lets kernels write
`eigen(mpts.s.ϵn[p])` directly on a stored strain object.
"""
@inline LinearAlgebra.eigen(strain::AbstractStrain{S,T,L}) where {S,T,L} = eigen(Symmetric(get_tensor(strain)))

"""
    _trial_elastic_strain(ΔFᵢⱼ, strain::LogarithmicStrain) -> LogarithmicStrain

Push the stored logarithmic strain forward through the incremental deformation
gradient `ΔFᵢⱼ` and return the trial elastic logarithmic strain.
"""
@inline function _trial_elastic_strain(ΔFᵢⱼ::SMatrix{S,S,T,L}, strain::LogarithmicStrain{S,T,L}) where {S,T,L}
    # trial left Cauchy-Green tensor
    λ, n = eigen(strain)
    bᵢⱼ  = ΔFᵢⱼ * (n * diagm(exp.(T(2.0) * λ)) * n') * ΔFᵢⱼ'
    # logarithmic strain tensor from the eigen-decomposition of bᵢⱼ
    λ, n = eigen(Symmetric(bᵢⱼ))
    ϵᵢⱼ  = T(0.5) * (n * diagm(log.(λ)) * n')
    vol  = tr(ϵᵢⱼ) / T(3.0)
    dev  = ϵᵢⱼ - vol * SMatrix{S,S,T,L}(I)
    return LogarithmicStrain(vol, SMatrix{S,S,T,L}(dev))
end

"""
    _infinitesimal_strain(ΔFᵢⱼ::SMatrix{S,S,T,L}) -> InfinitesimalStrain

Small-strain tensor `ϵ = ½(ΔF + ΔFᵀ) - I`, split into `vol = tr(ϵ)/3` and `dev = ϵ - vol·I`.
"""
@inline function InfinitesimalStrain(ϵᵢⱼ::SMatrix{S,S,T,L}) where {S,T,L}
    vol = tr(ϵᵢⱼ) / T(3.0)
    return InfinitesimalStrain(vol, SMatrix{S,S,T,L}(ϵᵢⱼ - vol * SMatrix{S,S,T,L}(I)))
end
@inline function _infinitesimal_strain(ΔFᵢⱼ::SMatrix{S,S,T,L}) where {S,T,L}
    return InfinitesimalStrain(T(0.5) .* (ΔFᵢⱼ + ΔFᵢⱼ') .- SMatrix{S,S,T,L}(I))
end

"""
    InfinitesimalStrain(ϵ::SVector{3,T}) / InfinitesimalStrain(ϵ::SVector{6,T})

Build directly from a full engineering-Voigt strain vector, same convention as
`LogarithmicStrain(ϵ::SVector)`.
"""
@inline function InfinitesimalStrain(ϵ::SVector{3,T}) where {T}
    vol = (ϵ[1] + ϵ[2]) / T(3.0)
    return InfinitesimalStrain(vol, SMatrix{2,2,T,4}(ϵ[1] - vol, ϵ[3]/T(2.0), ϵ[3]/T(2.0), ϵ[2] - vol))
end
@inline function InfinitesimalStrain(ϵ::SVector{6,T}) where {T}
    vol = (ϵ[1] + ϵ[2] + ϵ[3]) / T(3.0)
    return InfinitesimalStrain(vol, SMatrix{3,3,T,9}(
        ϵ[1] - vol, ϵ[6]/T(2.0), ϵ[5]/T(2.0),
        ϵ[6]/T(2.0), ϵ[2] - vol, ϵ[4]/T(2.0),
        ϵ[5]/T(2.0), ϵ[4]/T(2.0), ϵ[3] - vol,
    ))
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Stress tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    AbstractStress{S,T,L} <: AbstractTensor{T}

Second-order stress tensor stored as `p::T` (pressure, **positive in compression**) +
`dev::SMatrix{S,S,T,L}` (deviatoric part), so the full tensor is `dev - p·I`.
Declared in `types/abstract.jl`; concrete subtypes and methods live here.
"""
AbstractStress

"""
    get_tensor(stress::AbstractStress{S,T,L}) -> SMatrix{S,S,T,L}

Reassemble the full stress tensor `dev - p·I` (`p` positive in compression).
"""
@inline function get_tensor(stress::AbstractStress{S,T,L}) where {S,T,L}
    return stress.dev - stress.p * SMatrix{S,S,T,L}(I)
end

"""
    get_voigt(stress::AbstractStress) -> SVector{3,T} (2D) / SVector{6,T} (3D)

Voigt stress vector `(σxx,σyy,σxy)` / `(σxx,σyy,σzz,σyz,σxz,σxy)` — no shear doubling
(unlike strain).
"""
@inline function get_voigt(stress::AbstractStress{2,T,L}) where {T,L}
    return SVector{3,T}(
        stress.dev[1,1] - stress.p,
        stress.dev[2,2] - stress.p,
        stress.dev[1,2],
    )
end
@inline function get_voigt(stress::AbstractStress{3,T,L}) where {T,L}
    return SVector{6,T}(
        stress.dev[1,1] - stress.p,
        stress.dev[2,2] - stress.p,
        stress.dev[3,3] - stress.p,
        stress.dev[2,3],
        stress.dev[1,3],
        stress.dev[1,2],
    )
end

"""
    CauchyStress{S,T,L} <: AbstractStress{S,T,L}

Cauchy (true) stress. Stored by `mpts.s.σᵢ`/`mpts.s.σn`.
"""
struct CauchyStress{S,T,L} <: AbstractStress{S,T,L}
    p  ::T
    dev::SMatrix{S,S,T,L}
    function CauchyStress{S,T,L}(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
    function CauchyStress(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
end
@adapt_struct CauchyStress

# dev given as a Voigt deviator only (no pressure part)
@inline function CauchyStress(p::T, dev::SVector{3,T}) where {T}
    return CauchyStress(p, SMatrix{2,2,T,4}(
        dev[1], dev[3],
        dev[3], dev[2],
    ))
end
# 3D: SMatrix fills column-major — dev = [dev1 dev4 dev5; dev4 dev2 dev6; dev5 dev6 dev3].
@inline function CauchyStress(p::T, dev::SVector{6,T}) where {T}
    return CauchyStress(p, SMatrix{3,3,T,9}(
        dev[1], dev[4], dev[5],
        dev[4], dev[2], dev[6],
        dev[5], dev[6], dev[3],
    ))
end
# full Voigt stress vector (not a deviator) — 2D, p = -tr(σ)/2
@inline function CauchyStress(σ::SVector{3,T}) where {T}
    p = -(σ[1] + σ[2]) / T(2.0)
    return CauchyStress(p, SMatrix{2,2,T,4}(
        σ[1] + p, σ[3],
        σ[3],     σ[2] + p,
    ))
end
# full Voigt stress vector (not a deviator) — 3D, p = -tr(σ)/3
@inline function CauchyStress(σ::SVector{6,T}) where {T}
    p = -(σ[1] + σ[2] + σ[3]) / T(3.0)
    return CauchyStress(p, SMatrix{3,3,T,9}(
        σ[1] + p, σ[6],     σ[5],
        σ[6],     σ[2] + p, σ[4],
        σ[5],     σ[4],     σ[3] + p,
    ))
end

"""
    KirchhoffStress{S,T,L} <: AbstractStress{S,T,L}

Kirchhoff stress `τ = J·σ`. Stored by `mpts.s.τᵢ`; `transform` divides by `J` to get
`mpts.s.σᵢ`.
"""
struct KirchhoffStress{S,T,L} <: AbstractStress{S,T,L}
    p  ::T
    dev::SMatrix{S,S,T,L}
    function KirchhoffStress{S,T,L}(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
    function KirchhoffStress(p::T, dev::SMatrix{S,S,T,L}) where {S,T,L}
        return new{S,T,L}(p, dev)
    end
end
@adapt_struct KirchhoffStress

# dev given as a Voigt deviator only (no pressure part)
@inline function KirchhoffStress(p::T, dev::SVector{3,T}) where {T}
    return KirchhoffStress(p, SMatrix{2,2,T,4}(
        dev[1], dev[3],
        dev[3], dev[2],
    ))
end
# 3D: SMatrix fills column-major — dev = [dev1 dev4 dev5; dev4 dev2 dev6; dev5 dev6 dev3].
@inline function KirchhoffStress(p::T, dev::SVector{6,T}) where {T}
    return KirchhoffStress(p, SMatrix{3,3,T,9}(
        dev[1], dev[4], dev[5],
        dev[4], dev[2], dev[6],
        dev[5], dev[6], dev[3],
    ))
end
# full Voigt stress vector (not a deviator) — 2D, p = -tr(τ)/2
@inline function KirchhoffStress(τ::SVector{3,T}) where {T}
    p = -(τ[1] + τ[2]) / T(2.0)
    return KirchhoffStress(p, SMatrix{2,2,T,4}(
        τ[1] + p, τ[3],
        τ[3],     τ[2] + p,
    ))
end
# full Voigt stress vector (not a deviator) — 3D, p = -tr(τ)/3
@inline function KirchhoffStress(τ::SVector{6,T}) where {T}
    p = -(τ[1] + τ[2] + τ[3]) / T(3.0)
    return KirchhoffStress(p, SMatrix{3,3,T,9}(
        τ[1] + p, τ[6],     τ[5],
        τ[6],     τ[2] + p, τ[4],
        τ[5],     τ[4],     τ[3] + p,
    ))
end

"""
    _trial_elastic_stress(strain::LogarithmicStrain, Kc, Gc) -> KirchhoffStress

Isotropic linear-elastic trial Kirchhoff stress: `p = -3·Kc·ϵvol` (positive in
compression), `dev = 2·Gc·ϵdev`.
"""
@inline function _trial_elastic_stress(strain::LogarithmicStrain{S,T,L}, Kc::T, Gc::T) where {S,T,L}
    P   = -T(3.0) * Kc * strain.vol
    dev =  T(2.0) * Gc * strain.dev
    return KirchhoffStress(P, dev)
end

Base.zero(::Type{CauchyStress{S,T,L}})    where {S,T,L} = CauchyStress(   zero(T), zero(SMatrix{S,S,T,L}))
Base.zero(::Type{KirchhoffStress{S,T,L}}) where {S,T,L} = KirchhoffStress(zero(T), zero(SMatrix{S,S,T,L}))

"""
    CauchyStress(τ::KirchhoffStress, J⁻¹) -> CauchyStress
    KirchhoffStress(σ::CauchyStress) -> KirchhoffStress

Cross-type conversions (`τ ↔ σ`), used by `transform` and the u-P `dynamic_relaxation`
path. Defined outside both structs since `CauchyStress`'s version needs `KirchhoffStress`,
which isn't defined until later in the file.
"""
@inline function CauchyStress(τ::KirchhoffStress{S,T,L}, J⁻¹::T) where {S,T,L}
    return CauchyStress(τ.p * J⁻¹, τ.dev * J⁻¹)
end
@inline function KirchhoffStress(σ::CauchyStress{S,T,L}) where {S,T,L}
    return KirchhoffStress(σ.p, σ.dev)
end

"""
    get_J2(stress::AbstractStress) -> T

Second deviatoric stress invariant `J₂ = ½ sᵢⱼsᵢⱼ`, re-derived from `get_voigt` (see
`AbstractTensor`) rather than trusting the stored `dev`.
"""
@inline function get_J2(stress::AbstractStress{2,T,L}) where {T,L}
    σ = get_voigt(stress)
    P = (σ[1] + σ[2]) / T(2.0)
    return T(0.5) * ((σ[1] - P)^2 + (σ[2] - P)^2) + σ[3]^2
end
@inline function get_J2(stress::AbstractStress{3,T,L}) where {T,L}
    σ = get_voigt(stress)
    P = (σ[1] + σ[2] + σ[3]) / T(3.0)
    return T(0.5) * ((σ[1] - P)^2 + (σ[2] - P)^2 + (σ[3] - P)^2) + σ[4]^2 + σ[5]^2 + σ[6]^2
end

"""
    get_τII(stress::AbstractStress) -> T

Second invariant `τII = √J₂` of the deviatoric stress.
"""
@inline function get_τII(stress::AbstractStress)
    return sqrt(get_J2(stress))
end
