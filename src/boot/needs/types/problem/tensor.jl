# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Strain and stress tensors: supertypes, concrete types, readers and constructors. Both are stored
# as a scalar (volumetric strain / pressure) plus a trace-free 3D Voigt deviator.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export AbstractStrain, InfinitesimalStrain, LogarithmicStrain
export AbstractStress, CauchyStress, KirchhoffStress
export get_tensor, get_J2, get_τII

# ------------------------------------------------ strain

"""
    AbstractStrain{T2}

Second-order strain tensor stored on `PointSolidPhase` (`ϵᵢ`/`ϵn`) as its volumetric part
`vol::T2 = tr(ϵ)` and its trace-free deviator `dev::SVector{6,T2}`, a 3D Voigt vector
`(xx, yy, zz, γyz, γxz, γxy)` in every dimension: a 2D problem is plane strain (`γyz = γxz = 0`,
`zz` carried), a 1D problem is uniaxial strain. Shears are stored in the **engineering** convention
(`γ = 2ϵ`), which is what the elastic stiffness `Del` acts on. `AbstractStress` (below) is stored the same
way (pressure and deviator, without shear doubling).

Kinematics works on tensors, so the two meet in one place each way: `get_tensor(ϵ) = dev + vol/3·I`
(the symmetric 3×3, shears halved; `get_tensor(ϵ, Val(D))` is its in-plane `D×D` block), and the
3×3 constructor, which splits a tensor into `vol = tr` and `dev`. Code that builds `(vol, dev)`
directly (the return maps) passes a trace-free `dev`.
"""
abstract type AbstractStrain{T2} end

@inline function get_tensor(t::AbstractStrain{T2}) where {T2}
    d, m, h = t.dev, t.vol / T2(3.0), T2(0.5)
    return SMatrix{3,3,T2,9}(d[1] + m, h*d[6], h*d[5], h*d[6], d[2] + m, h*d[4], h*d[5], h*d[4], d[3] + m)
end
@inline get_tensor(t::AbstractStrain, ::Val{D}) where {D} = get_tensor(t)[SOneTo(D), SOneTo(D)]


"""
    LinearAlgebra.eigen(strain::AbstractStrain)

Eigen-decomposition of the (symmetric) 3×3 strain tensor.
"""
@inline LinearAlgebra.eigen(t::AbstractStrain) = eigen(Symmetric(get_tensor(t)))

"""
    InfinitesimalStrain{T2} <: AbstractStrain{T2}

Small-strain tensor `ϵ = ½(ΔF + ΔFᵀ) - I`, used under `material.elastic == "hypoelastic"`.
Built from `(vol, dev)` or from the 3×3 tensor.
"""
struct InfinitesimalStrain{T2} <: AbstractStrain{T2}
    vol::T2
    dev::SVector{6,T2}
end
@adapt_struct InfinitesimalStrain

"""
    LogarithmicStrain{T2} <: AbstractStrain{T2}

Logarithmic (Hencky) elastic strain tensor, used under `material.elastic ∈
{"hencky","improved_hencky"}`. Built by `_trial_elastic_strain` (from the 3×3 tensor), or from
`(vol, dev)` by the return maps (`_elastic_strain`); both in `elasticity.jl`.
"""
struct LogarithmicStrain{T2} <: AbstractStrain{T2}
    vol::T2
    dev::SVector{6,T2}
end
@adapt_struct LogarithmicStrain

# from a symmetric 3×3 strain tensor, for either strain type: `vol = tr`, `dev = M − vol/3·I`
# (shears doubled to engineering)
@inline function (::Type{S})(M::SMatrix{3,3,T2,9}) where {S<:AbstractStrain,T2}
    vol = M[1,1] + M[2,2] + M[3,3]
    sph = vol / T2(3.0)
    return S(vol, SVector{6,T2}(M[1,1] - sph, M[2,2] - sph, M[3,3] - sph, T2(2.0)*M[2,3], T2(2.0)*M[1,3], T2(2.0)*M[1,2]))
end

Base.zero(::Type{InfinitesimalStrain{T2}}) where {T2} = InfinitesimalStrain(zero(T2), zero(SVector{6,T2}))
Base.zero(::Type{LogarithmicStrain{T2}})   where {T2} = LogarithmicStrain(zero(T2), zero(SVector{6,T2}))

# ------------------------------------------------ stress

"""
    AbstractStress{T2}

Second-order stress tensor stored on `PointSolidPhase` (`σᵢ`/`σn`/`τᵢ`) as its pressure
`p::T2 = -tr/3` (positive in compression) and its trace-free deviator `dev::SVector{6,T2}`, a 3D
Voigt vector `(xx, yy, zz, yz, xz, xy)` in every dimension, like `AbstractStrain` (above) but
without shear doubling: in 2D plane strain `σzz` is carried.

As for strain, tensors meet the split in one place each way: `get_tensor(σ) = dev - p·I` (symmetric
3×3; `get_tensor(σ, Val(D))` is its in-plane `D×D` block, for force assembly) and the 3×3
constructor. Other readers: `get_J2`/`get_τII`.
Code that builds `(p, dev)` directly (the elastic laws, the return maps) passes a trace-free `dev`.
"""
abstract type AbstractStress{T2} end

@inline function get_tensor(t::AbstractStress{T2}) where {T2}
    d, p = t.dev, t.p
    return SMatrix{3,3,T2,9}(d[1] - p, d[6], d[5], d[6], d[2] - p, d[4], d[5], d[4], d[3] - p)
end
@inline get_tensor(t::AbstractStress, ::Val{D}) where {D} = get_tensor(t)[SOneTo(D), SOneTo(D)]

"""
    get_J2(stress::AbstractStress) -> T2

Second deviatoric stress invariant `J₂ = ½ sᵢⱼsᵢⱼ` of the full 3D stress.
"""
@inline function get_J2(t::AbstractStress{T2}) where {T2}
    s = t.dev
    return T2(0.5) * (s[1]^2 + s[2]^2 + s[3]^2) + s[4]^2 + s[5]^2 + s[6]^2
end

"""
    get_τII(stress::AbstractStress) -> T2

Second invariant `τII = √J₂` of the deviatoric stress.
"""
@inline get_τII(t::AbstractStress) = sqrt(get_J2(t))

"""
    CauchyStress{T2} <: AbstractStress{T2}

Cauchy (true) stress. Stored by `mpts.s.σᵢ`/`mpts.s.σn`.
"""
struct CauchyStress{T2} <: AbstractStress{T2}
    p  ::T2
    dev::SVector{6,T2}
end
@adapt_struct CauchyStress

"""
    KirchhoffStress{T2} <: AbstractStress{T2}

Kirchhoff stress `τ = J·σ`. Stored by `mpts.s.τᵢ`; `transform` divides by `J` to get
`mpts.s.σᵢ`.
"""
struct KirchhoffStress{T2} <: AbstractStress{T2}
    p  ::T2
    dev::SVector{6,T2}
end
@adapt_struct KirchhoffStress

# from a symmetric 3×3 stress tensor: `p = -tr/3`, `dev = M + p·I`
@inline function CauchyStress(M::SMatrix{3,3,T2,9}) where {T2}
    p = -(M[1,1] + M[2,2] + M[3,3]) / T2(3.0)
    return CauchyStress(p, SVector{6,T2}(M[1,1] + p, M[2,2] + p, M[3,3] + p, M[2,3], M[1,3], M[1,2]))
end
"""
    CauchyStress(τ::KirchhoffStress, J⁻¹) -> CauchyStress

`σ = τ/J`, used by `transform` and the dynamic-relaxation `update`.
"""
@inline CauchyStress(τ::KirchhoffStress{T2}, J⁻¹::T2) where {T2} = CauchyStress(τ.p * J⁻¹, τ.dev * J⁻¹)

Base.zero(::Type{CauchyStress{T2}})    where {T2} = CauchyStress(zero(T2), zero(SVector{6,T2}))
Base.zero(::Type{KirchhoffStress{T2}}) where {T2} = KirchhoffStress(zero(T2), zero(SVector{6,T2}))

# ------------------------------------------------ both

# new in-plane D×D block `M`, out-of-plane components carried from `t`, for any stress or strain type
# (plane strain in 2D, uniaxial strain in 1D; in 3D `M` is the whole tensor): used by the Hencky
# push-forward, the hypoelastic strain increment (carrying zeros) and the Jaumann rotation
@inline function (::Type{S})(M::SMatrix{D,D,T2}, t::S) where {S<:Union{AbstractStrain,AbstractStress},D,T2}
    E = get_tensor(t)
    return S(SMatrix{3,3,T2,9}(ntuple(k -> (i = (k - 1) % 3 + 1; j = (k - 1) ÷ 3 + 1; i ≤ D && j ≤ D ? M[i,j] : E[i,j]), Val(9))))
end
