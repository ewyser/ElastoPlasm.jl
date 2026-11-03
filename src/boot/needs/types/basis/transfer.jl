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
    ApicTransfer{T2,D,L,VM} <: AbstractTransfer

Affine PIC (APIC) transfer scheme (`basis.trsfr = "apic"`). Unlike `StdTransfer` and `TpicTransfer`, it carries a
`D×D` matrix pair per point, used by APIC alone (`shpfun.jl`'s `Dij_nd`, `mapsto/fun/Bij.jl`'s `Bij`).

# Fields
- `Bᵢⱼ::VM`: the affine velocity matrix, an `SMatrix{D,D,T2,L}` per point.
- `Dᵢⱼ::VM`: the inertia-like matrix that normalises it.

`T2` is the float type, `D` the dimension, `L = D*D` (Julia cannot compute `D*D` in a field type, so it is a
parameter, computed by the constructor), `VM` the vectors' type (host or device).

# Constructors
    ApicTransfer{T2,D}(nmp)
    ApicTransfer{T2,D,L}(Bᵢⱼ, Dᵢⱼ)

The first, zero matrices for `nmp` points (`get_transfer`).
"""
struct ApicTransfer{T2,D,L,VM} <: AbstractTransfer
    Bᵢⱼ::VM
    Dᵢⱼ::VM
end
function ApicTransfer{T2,D}(nmp::Integer) where {T2,D}
    L = D*D
    return ApicTransfer{T2,D,L}([zero(SMatrix{D,D,T2,L}) for _ in 1:nmp], [zero(SMatrix{D,D,T2,L}) for _ in 1:nmp])
end
function ApicTransfer{T2,D,L}(Bᵢⱼ::VM, Dᵢⱼ::VM) where {T2,D,L,VM}
    return ApicTransfer{T2,D,L,VM}(Bᵢⱼ, Dᵢⱼ)
end
function Adapt.adapt_structure(to, tr::ApicTransfer{T2,D,L}) where {T2,D,L}
    return adapt_fields(to, tr, ApicTransfer{T2,D,L})
end

"""
    get_transfer(which::String, ::Type{T2}, D::Integer, nmp::Integer) -> AbstractTransfer

Map a solver's `basis.trsfr` string to its concrete `AbstractTransfer` instance.
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
