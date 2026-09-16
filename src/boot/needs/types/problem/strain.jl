# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Strain tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractStrain, InfinitesimalStrain, LogarithmicStrain
export get_tensor, get_voigt

"""
    AbstractStrain{S,T,L}

Second-order strain tensor stored on `PointSolidPhase` (`ϵᵢⱼ`/`ϵn`) as a volumetric+
deviatoric split: `vol::T` (volumetric part) + `dev::SMatrix{S,S,T,L}` (deviatoric
part). `S` is the spatial dimension, `L == S*S`. `stress.jl`'s `AbstractStress`
follows the same `(p, dev)`-split convention for stresses.

`dev` is **not** required to be trace-free: `_trial_elastic_stress` (`stress.jl`) uses
`tr(ϵ)/3` even in 2D (plane-strain), while the Voigt constructors below use
`-tr(σ)/S` for stress. Always read a tensor via `get_voigt`/`get_tensor`, never
`.dev`/`.p` directly — `get_J2`/`get_τII` (`stress.jl`) re-derive a trace-free
deviator from `get_voigt` for exactly this reason.

`get_voigt` is one function dispatched onto `AbstractStrain`/`AbstractStress`: strain
uses the engineering-shear convention (`γxy = 2ϵxy`), stress does not.
"""
abstract type AbstractStrain{S,T,L} end

Base.getindex(strain::AbstractStrain{S,T,L}, i::Int, j::Int) where {S,T,L} =
    strain.dev[i,j] + (i == j ? strain.vol / T(3.0) : zero(T))

"""
    get_tensor(strain::AbstractStrain{S,T,L}) -> SMatrix{S,S,T,L}

Reassemble the full strain tensor `dev + (vol/3)·I`.
"""
@inline function get_tensor(strain::AbstractStrain{S,T,L}) where {S,T,L}
    return strain.dev + (strain.vol / T(3.0)) * SMatrix{S,S,T,L}(I)
end

"""
    get_voigt(strain::AbstractStrain) -> SVector{3,T} (2D) / SVector{6,T} (3D)

Voigt strain vector `(ϵxx,ϵyy,γxy)` / `(ϵxx,ϵyy,ϵzz,γyz,γxz,γxy)`, engineering-shear
convention (`γxy = 2ϵxy`) — required for `Del`-facing algebra (`σ = Del·ε_voigt`).
Round-trips with `LogarithmicStrain(ϵ::SVector)`/`InfinitesimalStrain(ϵ::SVector)`.
"""
@inline function get_voigt(strain::AbstractStrain{2,T,L}) where {T,L}
    return SVector{3,T}(
        strain.dev[1,1] + strain.vol / T(3.0),
        strain.dev[2,2] + strain.vol / T(3.0),
        T(2.0) * strain.dev[1,2],
    )
end
@inline function get_voigt(strain::AbstractStrain{3,T,L}) where {T,L}
    return SVector{6,T}(
        strain.dev[1,1] + strain.vol / T(3.0),
        strain.dev[2,2] + strain.vol / T(3.0),
        strain.dev[3,3] + strain.vol / T(3.0),
        T(2.0) * strain.dev[2,3],
        T(2.0) * strain.dev[1,3],
        T(2.0) * strain.dev[1,2],
    )
end

"""
    LinearAlgebra.eigen(strain::AbstractStrain)

Eigen-decomposition of the reassembled (symmetric) strain tensor. Lets kernels write
`eigen(mpts.s.ϵn[p])` directly on a stored strain object.
"""
@inline LinearAlgebra.eigen(strain::AbstractStrain{S,T,L}) where {S,T,L} = eigen(Symmetric(get_tensor(strain)))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# InfinitesimalStrain
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    _infinitesimal_strain(ΔFᵢⱼ::SMatrix{S,S,T,L}) -> InfinitesimalStrain

Small-strain tensor `ϵ = ½(ΔF + ΔFᵀ) - I`, split into `vol = tr(ϵ)/3` and `dev = ϵ - vol·I`.
"""
@inline function InfinitesimalStrain(ϵᵢⱼ::SMatrix{S,S,T,L}) where {S,T,L}
    vol = tr(ϵᵢⱼ)
    return InfinitesimalStrain(vol, SMatrix{S,S,T,L}(ϵᵢⱼ - vol / T(3.0) * SMatrix{S,S,T,L}(I)))
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
    vol = (ϵ[1] + ϵ[2])
    return InfinitesimalStrain(vol, SMatrix{2,2,T,4}(
        ϵ[1] - vol / T(3.0), ϵ[3]/T(2.0)        ,
        ϵ[3]/T(2.0)        , ϵ[2] - vol / T(3.0),
    ))
end
@inline function InfinitesimalStrain(ϵ::SVector{6,T}) where {T}
    vol = (ϵ[1] + ϵ[2] + ϵ[3])
    return InfinitesimalStrain(vol, SMatrix{3,3,T,9}(
        ϵ[1] - vol / T(3.0), ϵ[6]/T(2.0)        , ϵ[5]/T(2.0),
        ϵ[6]/T(2.0)        , ϵ[2] - vol / T(3.0), ϵ[4]/T(2.0),
        ϵ[5]/T(2.0)        , ϵ[4]/T(2.0)        , ϵ[3] - vol / T(3.0),
    ))
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LogarithmicStrain
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    vol = tr(ϵᵢⱼ)
    return LogarithmicStrain(vol, SMatrix{S,S,T,L}(ϵᵢⱼ - vol / T(3.0) * SMatrix{S,S,T,L}(I)))
end

"""
    LogarithmicStrain(ϵ::SVector{3,T}) / LogarithmicStrain(ϵ::SVector{6,T})

Build directly from a full engineering-Voigt strain vector. `vol = tr(ϵ)/3` (same
plane-strain convention as the `SMatrix` constructor above — not stress's `tr/S`).
Exact inverse of `get_voigt`.
"""
@inline function LogarithmicStrain(ϵ::SVector{3,T}) where {T}
    vol = (ϵ[1] + ϵ[2])
    return LogarithmicStrain(vol, SMatrix{2,2,T,4}(ϵ[1] - vol / T(3.0), ϵ[3]/T(2.0), ϵ[3]/T(2.0), ϵ[2] - vol / T(3.0)))
end
@inline function LogarithmicStrain(ϵ::SVector{6,T}) where {T}
    vol = (ϵ[1] + ϵ[2] + ϵ[3])
    return LogarithmicStrain(vol, SMatrix{3,3,T,9}(
        ϵ[1] - vol / T(3.0), ϵ[6]/T(2.0)        , ϵ[5]/T(2.0)        ,
        ϵ[6]/T(2.0)        , ϵ[2] - vol / T(3.0), ϵ[4]/T(2.0)        ,
        ϵ[5]/T(2.0)        , ϵ[4]/T(2.0)        , ϵ[3] - vol / T(3.0),
    ))
end

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
    vol  = tr(ϵᵢⱼ)
    dev  = ϵᵢⱼ - vol / T(3.0) * SMatrix{S,S,T,L}(I)
    return LogarithmicStrain(vol, SMatrix{S,S,T,L}(dev))
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Shared
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Base.zero(::Type{InfinitesimalStrain{S,T,L}}) where {S,T,L} = InfinitesimalStrain(zero(T), zero(SMatrix{S,S,T,L}))
Base.zero(::Type{LogarithmicStrain{S,T,L}})   where {S,T,L} = LogarithmicStrain(  zero(T), zero(SMatrix{S,S,T,L}))
