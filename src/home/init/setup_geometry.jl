"""
    Geometry(L::Vector{T2}, nel::Vector{T1}, instr::Instruction{T1,T2,D}; x₀::Vector=...) where {T1,T2,D}

Constructor for the `Geometry` type. Computes mesh geometry parameters and returns a `Geometry` struct containing mesh and domain information, based on the number of elements, domain size, and basis type.

# Arguments
- `nel::Vector{T1}`: Number of elements in each spatial direction.
- `L::Vector{T2}`: Length of the domain in each spatial direction.
- `instr`: Instruction object or dictionary containing basis information.

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
geom = get_geom([10, 10], [1.0, 1.0], Dict(:basis => Dict(:which => "bsmpm", :ghost => 1)))
@show geom.dim, geom.h, geom.nel, geom.nno, geom.L, geom.nn, geom.xB
```
"""

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0,]) where {T1<:Integer,T2<:Real,D<:OneDimension, S<:AbstractSolver{T1,T2,D}}    
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel
    
    # Add ghost element(s) (?) and modify nel vector accordingly
    ghost = solver.basis.ghost
    nel   = [nel[1]+ghost, (nel[1]+ghost),]
    xB    = [x₀[1]-ghost*h[1] L[1]+ghost*h[1]]

    # Create nno vector
    nno = [nel[1]+1, nel[1]+1,]
    
    # Define shape function compact support based on basis type
    if solver.basis.which == "smpm"
        neighbour  = neighbs((T1(0):T1(1),))
    else
        neighbour  = neighbs((T1(-1):T1(2),))
    end

    return Geometry{T1,T2,D,typeof(neighbour)}(
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

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0, 0.0,]) where {T1<:Integer,T2<:Real,D<:TwoDimension, S<:AbstractSolver{T1,T2,D}}    
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    # Add ghost element(s) and modify nel vector accordingly
    ghost = solver.basis.ghost
    nel   = [nel[1]+ghost, nel[2]+ghost, (nel[1]+ghost) * (nel[2]+ghost),]
    xB    = vcat([x₀[1]-ghost*h[1] L[1]+ghost*h[1]], [x₀[2]-ghost*h[2] L[2]+ghost*h[2]])    

    # Create nno vector
    nno = [nel[1]+1, nel[2]+1, (nel[1]+1) * (nel[2]+1),]
    
    # Define shape function compact support based on basis type
    if solver.basis.which == "smpm"
        neighbour  = neighbs((T1(0):T1(1), T1(0):T1(1)))
    else
        neighbour  = neighbs((T1(-1):T1(2), T1(-1):T1(2)))
    end

    return Geometry{T1,T2,D,typeof(neighbour)}(
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

function Geometry(L::Vector, nel::Vector, solver::S; x₀::Vector=[0.0, 0.0, 0.0]) where {T1<:Integer,T2<:Real,D<:ThreeDimension, S<:AbstractSolver{T1,T2,D}}
    # Calculate problem dimensionality & node spacing
    dim,h = length(L),L ./ nel

    # Add ghost element(s) (?) and modify nel vector accordingly
    ghost = instr.basis.ghost
    nel   = [nel[1]+ghost, nel[2]+ghost, nel[3]+ghost, (nel[1]+ghost) * (nel[2]+ghost) * (nel[3]+ghost),]
    xB    = vcat([x₀[1]-ghost*h[1] L[1]+ghost*h[1]], [x₀[2]-ghost*h[2] L[2]+ghost*h[2]], [x₀[3]-ghost*h[3] L[3]+ghost*h[3]])  

    # Create nno vector
    nno = [nel[1]+1, nel[2]+1, nel[3]+1, (nel[1]+1) * (nel[2]+1) * (nel[3]+1),]
    
    # Define shape function compact support based on basis type
    if instr.basis.which == "smpm"
        neighbour  = neighbs((T1(0):T1(1), T1(0):T1(1), T1(0):T1(1)))
    else
        neighbour  = neighbs((T1(-1):T1(2), T1(-1):T1(2), T1(-1):T1(2)))
    end

    return Geometry{T1,T2,D,typeof(neighbour)}(
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

function setup_geometry(L::Vector{T2}, nel::Vector{T1}, solver::S) where {T1,T2,D,S<:AbstractSolver}
    return Geometry(L,nel,solver)
end