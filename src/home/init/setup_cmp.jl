"""
    get_elastic_moduli(E, ν) -> Kc, Gc

Bulk and shear moduli from Young's modulus and Poisson's ratio. The elastic laws use them directly;
the Voigt stiffness `Del` can be built from them by `get_Del` (`elasticity.jl`).
"""
function get_elastic_moduli(E::T2,ν::T2) where {T2}
    return T2(E/(3.0*(1.0-2.0*ν))), T2(E/(2.0*(1.0+ν)))
end

"""
    setup_cmp(nmp, coh0, cohr, phi; E, ν, Hp, D, constitutive) -> Vector{<:AbstractConstitutiveModel}

The constitutive constants of each material point: `DruckerPrager` (`constitutive = "DP"`) or
`VonMises` (`"VM"`, which has no friction angle, so `phi` is ignored) from its cohesion `coh0`,
residual cohesion `cohr`, friction angle `phi`, softening modulus `Hp` and the bulk and shear moduli
of `E`, `ν` (`get_elastic_moduli`); `E`, `ν`, `Hp` are one value or one per point. Any other
`constitutive` is an error here, at setup. The elastic law (`material.elastic`) is not part of `cmp`:
it is `Point`'s `EL` type parameter, set by `build_solid_phase`.
"""
function setup_cmp(nmp::Integer,coh0::AbstractVector{T2},cohr::AbstractVector{T2},phi::AbstractVector{T2};
    E::Union{T2,AbstractVector{T2}},ν::Union{T2,AbstractVector{T2}},Hp::Union{T2,AbstractVector{T2}},D::Integer,constitutive::String,
    ) where {T2}
    each(x) = x isa AbstractVector ? x : fill(x, nmp)
    KG, H = get_elastic_moduli.(each(E), each(ν)), each(Hp)
    if constitutive == "DP"
        return [DruckerPrager{T2,D}(KG[p][2],KG[p][1],H[p],T2(coh0[p]),T2(cohr[p]),T2(phi[p])) for p ∈ 1:nmp]
    elseif constitutive == "VM"
        return [VonMises{T2,D}(KG[p][2],KG[p][1],H[p],T2(coh0[p]),T2(cohr[p])) for p ∈ 1:nmp]
    else
        throw(error("InvalidConstitutiveModel: $(constitutive) (expected \"DP\" or \"VM\")"))
    end
end