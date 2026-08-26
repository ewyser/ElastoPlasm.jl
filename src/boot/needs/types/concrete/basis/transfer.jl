# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete AbstractTransfer types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export StdTransfer, TpicTransfer, ApicTransfer, get_transfer

"""
    StdTransfer

P2G/G2P transfer scheme with no extra per-particle state — plain PIC/FLIP-blended
transfer. Carries no fields.
"""
struct StdTransfer <: AbstractTransfer end

"""
    TpicTransfer

Technique-corrected PIC (TPIC) transfer scheme. Carries no fields — the per-particle
velocity-gradient correction it uses (`mpts.s.∇vᵢⱼ`) already lives on `Point`.
"""
struct TpicTransfer <: AbstractTransfer end

"""
    ApicTransfer{T,D,L}

Affine PIC (APIC) transfer scheme. Unlike `StdTransfer`/`TpicTransfer`, APIC needs two
per-particle `D×D` accumulators (`Bᵢⱼ`, the affine velocity matrix, and `Dᵢⱼ`, the
inertia-like normalization matrix) — moved here from `Point` since they are 100%
APIC-exclusive (see `shpfun.jl`'s `Dij_nd` and `mapsto/fun/Bij.jl`'s `Bij`, both only
registered/called under `transfer.trsfr=="apic"`). Needs `nmp` at construction time to
size its vectors, unlike `StdTransfer`/`TpicTransfer`'s zero-arg constructors. `L=D*D`
is a separate type parameter (matching `DruckerPrager{T2,D,NSTR,L}`'s pattern) since
`SMatrix{D,D,T,D*D}` — computing the length inline from `D` — isn't valid inside a
struct's own field-type declaration (`D*D` on a bare `TypeVar`); `L` is only computed
once `D` is a concrete value, inside the inner constructor.
"""
struct ApicTransfer{T,D,L} <: AbstractTransfer
    Bᵢⱼ::Vector{SMatrix{D,D,T,L}}
    Dᵢⱼ::Vector{SMatrix{D,D,T,L}}
    #
    function ApicTransfer{T,D}(nmp::Integer) where {T,D}
        L = D*D
        return new{T,D,L}(
            [zero(SMatrix{D,D,T,L}) for _ in 1:nmp],
            [zero(SMatrix{D,D,T,L}) for _ in 1:nmp],
        )
    end
end

"""
    get_transfer(which::String, ::Type{T2}, D::Integer, nmp::Integer) -> AbstractTransfer

Map a solver's `transfer.trsfr` string to its concrete `AbstractTransfer` instance.
Mirrors `get_basis`, except `ApicTransfer` additionally needs `T2`/`D`/`nmp` to size its
per-particle vectors (which `get_basis`'s marker-only kinds never needed).
"""
function get_transfer(which::String, ::Type{T2}, D::Integer, nmp::Integer) where {T2}
    if which == "std"
        return StdTransfer()
    elseif which == "tpic"
        return TpicTransfer()
    elseif which == "apic"
        return ApicTransfer{T2,D}(nmp)
    else
        return error("Unsupported transfer scheme: $which")
    end
end
