# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Strain tensor supertype, concrete types and methods
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractStrain, InfinitesimalStrain, LogarithmicStrain
export get_tensor

"""
    AbstractStrain{T}

Second-order strain tensor stored on `PointSolidPhase` (`ϵᵢ`/`ϵn`) as its volumetric part
`vol::T = tr(ϵ)` and its trace-free deviator `dev::SVector{6,T}`, a 3D Voigt vector
`(xx, yy, zz, γyz, γxz, γxy)` in every dimension: a 2D problem is plane strain (`γyz = γxz = 0`,
`zz` carried), a 1D problem is uniaxial strain. Shears are stored in the **engineering** convention
(`γ = 2ϵ`), which is what the elastic stiffness `Del` acts on. `stress.jl`'s `AbstractStress` is
stored the same way (pressure and deviator, without shear doubling).

Kinematics works on tensors, so the two meet in one place each way: `get_tensor(ϵ) = dev + vol/3·I`
(the symmetric 3×3, shears halved; `get_tensor(ϵ, Val(D))` is its in-plane `D×D` block), and the
3×3 constructor, which splits a tensor into `vol = tr` and `dev`. Code that builds `(vol, dev)`
directly (the return maps) passes a trace-free `dev`.
"""
abstract type AbstractStrain{T} end

@inline function get_tensor(t::AbstractStrain{T}) where {T}
    d, m, h = t.dev, t.vol / T(3.0), T(0.5)
    return SMatrix{3,3,T,9}(d[1] + m, h*d[6], h*d[5], h*d[6], d[2] + m, h*d[4], h*d[5], h*d[4], d[3] + m)
end
@inline get_tensor(t::AbstractStrain, ::Val{D}) where {D} = get_tensor(t)[SOneTo(D), SOneTo(D)]


"""
    LinearAlgebra.eigen(strain::AbstractStrain)

Eigen-decomposition of the (symmetric) 3×3 strain tensor.
"""
@inline LinearAlgebra.eigen(t::AbstractStrain) = eigen(Symmetric(get_tensor(t)))

"""
    InfinitesimalStrain{T} <: AbstractStrain{T}

Small-strain tensor `ϵ = ½(ΔF + ΔFᵀ) - I`, used under `material.elastic == "hypoelastic"`.
Built from `(vol, dev)` or from the 3×3 tensor.
"""
struct InfinitesimalStrain{T} <: AbstractStrain{T}
    vol::T
    dev::SVector{6,T}
end
@adapt_struct InfinitesimalStrain

"""
    LogarithmicStrain{T} <: AbstractStrain{T}

Logarithmic (Hencky) elastic strain tensor, used under `material.elastic ∈
{"hencky","improved_hencky"}`. Built by `_trial_elastic_strain` (from the 3×3 tensor), or from
`(vol, dev)` by the return maps (`_elastic_strain`, `stress.jl`).
"""
struct LogarithmicStrain{T} <: AbstractStrain{T}
    vol::T
    dev::SVector{6,T}
end
@adapt_struct LogarithmicStrain

# from a symmetric 3×3 strain tensor, for either strain type: `vol = tr`, `dev = M − vol/3·I`
# (shears doubled to engineering)
@inline function (::Type{S})(M::SMatrix{3,3,T,9}) where {S<:AbstractStrain,T}
    vol = M[1,1] + M[2,2] + M[3,3]
    m   = vol / T(3.0)
    return S(vol, SVector{6,T}(M[1,1] - m, M[2,2] - m, M[3,3] - m, T(2.0)*M[2,3], T(2.0)*M[1,3], T(2.0)*M[1,2]))
end

"""
    _trial_elastic_strain(ΔFᵢⱼ::SMatrix{D,D}, ϵᵢ::LogarithmicStrain) -> LogarithmicStrain

Push the stored logarithmic strain forward through the incremental deformation gradient:
`bᵉ = ΔF·exp(2ϵ)·ΔFᵀ`, `ϵ = ½ ln bᵉ`. Only the in-plane `D×D` block is pushed forward: in plane
strain (2D) and uniaxial strain (1D) `ΔF` has no out-of-plane stretch or shear (`ΔF_zz = 1`), so
`bᵉ` stays block-diagonal and the out-of-plane components are unchanged (`b_zz = exp(2ϵzz)`, hence
`ϵzz` is kept as is); `LogarithmicStrain(ϵᵢⱼ, ϵᵢ)` builds the new strain from the pushed-forward
block, carrying the out-of-plane components (`stress.jl`).
"""
@inline function _trial_elastic_strain(ΔFᵢⱼ::SMatrix{D,D,T}, ϵᵢ::LogarithmicStrain{T}) where {D,T}
    # get the in-plane block of the previously-converged logarithmic strain
    ϵᵢⱼ  = get_tensor(ϵᵢ, Val(D))
    λ, n = eigen(Symmetric(ϵᵢⱼ))
    # calculate left Cauchy-Green strain
    bᵢⱼ  = ΔFᵢⱼ * (n * diagm(exp.(T(2.0) * λ)) * n') * ΔFᵢⱼ'
    # update logarithmic strain tensor from the eigen-decomposition of bᵢⱼ
    λ, n = eigen(Symmetric(bᵢⱼ))
    ϵᵢⱼ  = T(0.5) * (n * diagm(log.(λ)) * n')
    return LogarithmicStrain(ϵᵢⱼ, ϵᵢ)
end

Base.zero(::Type{InfinitesimalStrain{T}}) where {T} = InfinitesimalStrain(zero(T), zero(SVector{6,T}))
Base.zero(::Type{LogarithmicStrain{T}})   where {T} = LogarithmicStrain(zero(T), zero(SVector{6,T}))
