"""
    set_roller_dirichlet(nno, xn, xB)

Set roller Dirichlet boundary conditions for a mesh.

# Arguments
- `nno`: Number of nodes in each direction and total.
- `xn`: Matrix of nodal coordinates.
- `xB`: Boundary coordinates.

# Returns
- `bc::Matrix{Bool}`: Boolean matrix indicating roller boundary conditions.
"""
function set_dirichlet(bc,xn,lower_upper_lim,dim; node::Symbol=:roller)
    id = findall(x-> x ∈ lower_upper_lim,xn[dim,:])
    if node == :roller
        bc[dim,id].= true
    elseif node == :fixed
        bc[:  ,id].= true
    end
    return nothing
end

"""
    get_bc(xn::Matrix{T2}, instr::NamedTuple; ghosts::Vector{T2}=[T2(0.0)]) where {T2}

Compute the boundary condition matrix and boundary coordinates for a mesh, using boundary condition types specified in `instr.bcs.dirichlet`.

# Arguments
- `xn::Matrix{T2}`: Matrix of nodal coordinates (size: ndim × nno).
- `instr::NamedTuple`: Instruction/configuration named tuple, must contain `bcs.dirichlet` as a matrix of Symbols (e.g., `:roller`, `:fixed`).
- `ghosts::Vector{T2}`: Optional, size of ghost cells to add at boundaries (default: `[T2(0.0)]`).

# Returns
- `bc::Matrix{Bool}`: Boolean matrix (ndim × nno) indicating where boundary conditions are applied.
- `xB::Matrix{T2}`: Matrix (ndim × 2) of minimum and maximum boundary coordinates (with ghosts applied).

# Example
```julia
bc, xB = get_bc(xn, instr; ghosts=[0.0])
```
"""
function get_bc(xn::Matrix{T2},instr::NamedTuple; ghosts::Vector{T2}=[T2(0.0)]) where {T2}
    xB = hcat(
        minimum(xn,dims=2).+ghosts,
        maximum(xn,dims=2).-ghosts
    )
    ndim = size(xn,1)
    nno  = size(xn,2)
    bc   = zeros(Bool,ndim,nno[end])
    for dim ∈ 1:ndim
        for (k,limit) ∈ enumerate(xB[dim,:])
            set_dirichlet(bc,xn,limit,dim; node = instr.bcs.dirichlet[dim,k]) 
        end
    end
    return bc,xB
end