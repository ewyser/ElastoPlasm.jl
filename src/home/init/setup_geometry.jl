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
    - `nel`: Number of elements (possibly with ghost elements)
    - `nno`: Number of nodes (possibly with ghost nodes)
    - `L`: Domain size in each direction
    - `nn`: Neighbourhood information for shape functions
    - `xB`: Domain bounds (with ghost region)

# Example
```julia
geom = Geometry([1.0, 1.0], [10, 10], solver)
@show geom.dim, geom.h, geom.nel, geom.nno, geom.L, geom.nn, geom.xB
```
"""

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0,]) where {T1<:Integer,T2<:Real, S<:AbstractSolver{T1,T2,1}}
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    # Add ghost element(s) (?) and modify nel vector accordingly
    ghost = solver.basis.ghost
    nel   = [nel[1]+ghost, (nel[1]+ghost),]
    xB    = [x₀[1]-ghost*h[1] L[1]+ghost*h[1]]

    # Create nno vector
    nno = [nel[1]+1, nel[1]+1,]

    # Define shape function compact support based on basis type
    neighbour = Neighbs(get_basis(solver.basis.which), dim, T1)
    NN = neighbour.nn

    return Geometry{T1,T2,1,NN}(
        T1(dim),
        T1(ghost),
        T2.(h),
        T1.(nel),
        T1.(nno),
        T2.(L),
        neighbour,
        T2.(xB),
    )
end

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0, 0.0,]) where {T1<:Integer,T2<:Real, S<:AbstractSolver{T1,T2,2}}
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    # Add ghost element(s) and modify nel vector accordingly
    ghost = solver.basis.ghost
    nel   = [nel[1]+ghost, nel[2]+ghost, (nel[1]+ghost) * (nel[2]+ghost),]
    xB    = vcat([x₀[1]-ghost*h[1] L[1]+ghost*h[1]], [x₀[2]-ghost*h[2] L[2]+ghost*h[2]])

    # Create nno vector
    nno = [nel[1]+1, nel[2]+1, (nel[1]+1) * (nel[2]+1),]

    # Define shape function compact support based on basis type
    neighbour = Neighbs(get_basis(solver.basis.which), dim, T1)
    NN = neighbour.nn

    return Geometry{T1,T2,2,NN}(
        T1(dim),
        T1(ghost),
        T2.(h),
        T1.(nel),
        T1.(nno),
        T2.(L),
        neighbour,
        T2.(xB),
    )
end

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0, 0.0, 0.0]) where {T1<:Integer,T2<:Real, S<:AbstractSolver{T1,T2,3}}
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    # Add ghost element(s) (?) and modify nel vector accordingly
    ghost = solver.basis.ghost
    nel   = [nel[1]+ghost, nel[2]+ghost, nel[3]+ghost, (nel[1]+ghost) * (nel[2]+ghost) * (nel[3]+ghost),]
    xB    = vcat([x₀[1]-ghost*h[1] L[1]+ghost*h[1]], [x₀[2]-ghost*h[2] L[2]+ghost*h[2]], [x₀[3]-ghost*h[3] L[3]+ghost*h[3]])

    # Create nno vector
    nno = [nel[1]+1, nel[2]+1, nel[3]+1, (nel[1]+1) * (nel[2]+1) * (nel[3]+1),]

    # Define shape function compact support based on basis type
    neighbour = Neighbs(get_basis(solver.basis.which), dim, T1)
    NN = neighbour.nn

    return Geometry{T1,T2,3,NN}(
        T1(dim),
        T1(ghost),
        T2.(h),
        T1.(nel),
        T1.(nno),
        T2.(L),
        neighbour,
        T2.(xB),
    )
end

function setup_geometry(L::Vector{T2}, nel::Vector{T1}, solver::S) where {T1,T2,S<:AbstractSolver}
    return Geometry(L,nel,solver)
end