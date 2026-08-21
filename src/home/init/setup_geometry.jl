"""
    Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=...) where {S<:AbstractSolver}

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

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0,]) where {T1<:Integer,T2<:Real, S<:AbstractSolver{T1,T2,1}}
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    nel = [nel[1], nel[1],]
    xB  = [x₀[1] L[1]]

    # Create nno vector
    nno = [nel[1]+1, nel[1]+1,]

    return Geometry{T1,T2,1}(
        T1(dim),
        T2.(h),
        T1.(nel),
        T1.(nno),
        T2.(L),
        T2.(xB),
    )
end

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0, 0.0,]) where {T1<:Integer,T2<:Real, S<:AbstractSolver{T1,T2,2}}
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    nel = [nel[1], nel[2], nel[1]*nel[2],]
    xB  = vcat([x₀[1] L[1]], [x₀[2] L[2]])

    # Create nno vector
    nno = [nel[1]+1, nel[2]+1, (nel[1]+1) * (nel[2]+1),]

    return Geometry{T1,T2,2}(
        T1(dim),
        T2.(h),
        T1.(nel),
        T1.(nno),
        T2.(L),
        T2.(xB),
    )
end

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0, 0.0, 0.0]) where {T1<:Integer,T2<:Real, S<:AbstractSolver{T1,T2,3}}
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    nel = [nel[1], nel[2], nel[3], nel[1]*nel[2]*nel[3],]
    xB  = vcat([x₀[1] L[1]], [x₀[2] L[2]], [x₀[3] L[3]])

    # Create nno vector
    nno = [nel[1]+1, nel[2]+1, nel[3]+1, (nel[1]+1) * (nel[2]+1) * (nel[3]+1),]

    return Geometry{T1,T2,3}(
        T1(dim),
        T2.(h),
        T1.(nel),
        T1.(nno),
        T2.(L),
        T2.(xB),
    )
end

function setup_geometry(L::Vector{T2}, nel::Vector{T1}, solver::S) where {T1,T2,S<:AbstractSolver}
    return Geometry(L,nel,solver)
end
