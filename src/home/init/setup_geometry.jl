"""
    Geometry(L::Vector, nel::Vector, solver::AbstractSolver{T1,T2,D}; x₀=zeros(length(L))) -> Geometry{T1,T2,D}

Constructor for the `Geometry` type. Computes mesh geometry parameters and returns a `Geometry` struct containing mesh and domain information, based on the number of elements, domain size, and basis type.

# Arguments
- `L::Vector`: Length of the domain in each spatial direction.
- `nel::Vector`: Number of elements in each spatial direction.
- `solver::S`: Solver instance (e.g. `ExplicitSolver`) containing basis information.

# Returns
- `Geometry` struct with fields:
    - `dim`: Number of spatial dimensions
    - `h`: Element size in each direction
    - `nel`: Number of elements
    - `nno`: Number of nodes
    - `L`: Domain size in each direction
    - `xB`: Domain bounds

# Example
```julia
geom = Geometry([1.0, 1.0], [10, 10], solver)
@show geom.dim, geom.h, geom.nel, geom.nno, geom.L, geom.xB
```
"""

function Geometry(L::Vector, nel::Vector, solver::AbstractSolver{T1,T2,D}; x₀::Vector=zeros(length(L))) where {T1,T2,D}
    # per-axis counts followed by their product (the total); in 1D this gives [n, n]
    nno = nel .+ 1
    return Geometry{T1,T2,D}(
        T1(length(L)),
        T2.(L ./ nel),
        T1.(vcat(nel, prod(nel))),
        T1.(vcat(nno, prod(nno))),
        T2.(L),
        T2.(hcat(x₀, L)),
    )
end

function setup_geometry(L::Vector{T2}, nel::Vector{T1}, solver::S) where {T1,T2,S<:AbstractSolver}
    return Geometry(L,nel,solver)
end
