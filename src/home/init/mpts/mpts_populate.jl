"""
    mpts_populate(props, mat, solver::S; ni=2) where {S<:AbstractSolver}

Initialize candidate material point fields and coordinates on a regular per-element grid, for use in MPM/ElastoPlasm simulations.

# Arguments
- `props`: `MeshProperties` object containing geometry and boundary information.
- `mat`: Material parameters (NamedTuple, see `setup_material_constants`), must include `:c0` and `:cr`.
- `solver::S`: Solver instance (e.g. `ExplicitSolver`); only its dimension is used.
- `ni`: Number of intervals per element (default: 2).

# Returns
- `NamedTuple`: Contains
    - `:ni`   — Number of intervals per element.
    - `:x`    — Matrix of material point coordinates (D×N).

# Example
```julia
fields = mpts_populate(props, mat, solver; ni=4)
@show fields.x
@show fields.c0
```
"""
function mpts_populate(props,mat,solver::S; ni = 2,) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    # ni points per element along each axis, at sub-cell centres
    axs = [collect(props.xB[d,1]+(0.5*props.h[d]/ni):props.h[d]/ni:props.xB[d,2]-(0.5*props.h[d]/ni)) for d ∈ 1:D]
    # point order: vertical (last) axis fastest, then x, then y
    pts = Iterators.product(axs[D], axs[1:D-1]...)
    x   = reshape(stack((t[2:end]..., t[1]) for t ∈ pts), D, :)
    return (; ni=ni, x=x)
end
