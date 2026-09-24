"""
    material_fields(mat, nmp; coh0=fill(mat[:c0],nmp), cr=mat[:cr], ϕ=mat[:ϕ0]) -> NamedTuple

Per-particle material fields for `nmp` particles: initial cohesion `coh0` (pass a vector, e.g. from
`property_field` for a random field), residual cohesion `cohr`, friction angle `phi`, and the
thermal fields `c`, `k`, `T`, uniform from `mat`. Problem-specific variations (e.g. `get_slump`'s
weaker base layer) are applied by the caller on the returned vectors.
"""
function material_fields(mat, nmp; coh0=fill(mat[:c0],nmp), cr=mat[:cr], ϕ=mat[:ϕ0])
    return (;
        coh0 = coh0,
        cohr = fill(cr, nmp),
        phi  = fill(ϕ, nmp),
        c    = fill(mat[:specific_heat_capacity], nmp),
        k    = fill(mat[:thermal_conductivity], nmp),
        T    = fill(mat[:initial_temperature], nmp),
    )
end
