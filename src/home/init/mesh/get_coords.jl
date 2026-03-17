"""
    get_coords(ndim::T1, L::Vector{T2}, h::Vector{T2}; ghosts::Vector{T2}=[T2(0.0)]) where {T1,T2}

Generate nodal coordinates and mesh topology for 1D, 2D, or 3D domains.

# Arguments
- `ndim::T1`: Number of spatial dimensions.
- `L::Vector{T2}`: Length of the domain in each spatial direction.
- `h::Vector{T2}`: Element size in each direction.
- `ghosts::Vector{T2}`: Optional, size of ghost cells to add at boundaries (default: `[T2(0.0)]`).

# Returns
- `x::Matrix{T2}`: Matrix of nodal coordinates (each column is a node).
- `nel::Vector{Int}`: Number of elements in each direction and total.
- `nno::Vector{Int}`: Number of nodes in each direction and total.

# Example
```julia
x, nel, nno = get_coords(2, [1.0, 1.0], [0.1, 0.1])
```
"""
function get_coords(geom)# where {T1,T2}
    ndim, h, xB = geom.dim, geom.h, geom.xB
    T1,T2 = eltype(ndim), eltype(h)


    if ndim == 1
        x   = collect(xB[1,1]:h[1]:xB[1,2])
        nno = [length(x),length(x)] 
        nel = [nno[1]-1 ,nno[1]-1 ]
        x   = vcat(vec(x))
    elseif ndim == 2
        xvec = collect(xB[1,1]:h[1]:xB[1,2])
        yvec = collect(xB[2,1]:h[2]:xB[2,2])
        nx, ny = length(xvec), length(yvec)
        nno = [nx, ny, nx*ny]
        nel = [nx-1, ny-1, (nx-1)*(ny-1)]

        # Create meshgrid in (nx, ny) order
        xmat = repeat(reshape(xvec, nx, 1), 1, ny)
        ymat = repeat(reshape(yvec, 1, ny), nx, 1)
        x = vcat(vec(xmat)', vec(ymat)')
    elseif ndim == 3
        xvec = collect(xB[1,1]:h[1]:xB[1,2])
        yvec = collect(xB[2,1]:h[2]:xB[2,2])
        zvec = collect(xB[3,1]:h[3]:xB[3,2])
        nx, ny, nz = length(xvec), length(yvec), length(zvec)
        nno = [nx, ny, nz, nx*ny*nz]
        nel = [nx-1, ny-1, nz-1, (nx-1)*(ny-1)*(nz-1)]

        # Create meshgrid in (x, y, z) order
        xmat = repeat(reshape(xvec, nx, 1, 1), 1, ny, nz)
        ymat = repeat(reshape(yvec, 1, ny, 1), nx, 1, nz)
        zmat = repeat(reshape(zvec, 1, 1, nz), nx, ny, 1)
        x = vcat(vec(xmat)', vec(ymat)', vec(zmat)')
    end
    return T2.(x),T1.(nel),T1.(nno)
end