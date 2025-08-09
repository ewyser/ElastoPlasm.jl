"""
    get_geom(nel::Vector{T1}, L::Vector{T2}, instr) where {T1,T2}

Compute mesh geometry parameters based on the number of elements, domain size, and basis type.

# Arguments
- `nel::Vector{T1}`: Number of elements in each spatial direction.
- `L::Vector{T2}`: Length of the domain in each spatial direction.
- `instr`: Instruction dictionary containing basis information.

# Returns
- `ndim::T1`: Number of spatial dimensions.
- `nn::T1`: Number of nodes.
- `h::Vector{T2}`: Element size in each direction.

# Example
```julia
ndim, nn, h = get_geom([10, 10], [1.0, 1.0], Dict(:basis => Dict(:which => "bsmpm")))
```
"""
function get_geom(nel::Vector{T1},L::Vector{T2},instr) where {T1,T2}
    if instr[:basis][:which] == "bsmpm"
        ndim,nn,h = length(L),4^length(L),min.(L./nel,L./4)
    else
        ndim,nn,h = length(L),4^length(L),(L./nel)
    end
    return T1(ndim),T1(nn),T2.(h)
end

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
function get_coords(ndim::T1,L::Vector{T2},h::Vector{T2}; ghosts::Vector{T2}=[T2(0.0)]) where {T1,T2}
    if ndim == 1
        x0  = [0.0-ghosts[1],L[1]+ghosts[1]]
        x   = collect(first(x0):h[1]:last(x0))
        nno = [length(x),length(x)] 
        nel = [nno[1]-1 ,nno[1]-1 ]
        x   = vcat(vec(x))
    elseif ndim == 2
        x0  = [0.0-ghosts[1],L[1]+ghosts[1]]
        z0  = [0.0-ghosts[2],L[2]+ghosts[2]]
        x,z = collect(first(x0):h[1]:last(x0)),collect(first(z0):h[2]:last(z0))
        nno = [length(x),length(z),length(x)*length(z)  ] 
        nel = [nno[1]-1 ,nno[2]-1 ,(nno[1]-1)*(nno[2]-1)]

        x   = repeat((reshape(x,1     ,nno[1])),nno[2],1     )
        z   = repeat((reshape(z,nno[2],1     )),     1,nno[1])
        x   = vcat(vec(x)',vec(z)')
    elseif ndim == 3
        x0  = [0.0-ghosts[1],L[1]+ghosts[1]]
        y0  = [0.0-ghosts[2],L[2]+ghosts[2]]
        z0  = [0.0-ghosts[3],L[3]+2.0*h[3]+ghosts[3]]
        x   = collect(first(x0):h[1]:last(x0)) 
        y   = collect(first(y0):h[2]:last(y0)) 
        z   = collect(first(z0):h[3]:last(z0))   
        nno = [length(x),length(y),length(z),length(x)*length(y)*length(z)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]

        x   = repeat((reshape(x,1     ,nno[1],1     )),nno[3],1     ,nno[2])
        y   = repeat((reshape(y,1     ,1     ,nno[2])),nno[3],nno[1],1     )
        z   = repeat((reshape(z,nno[3],1     ,1     )),1     ,nno[1],nno[2]) 
        x   = vcat(vec(x)',vec(y)',vec(z)')
    end
    return T2.(x),T1.(nel),T1.(nno)
end