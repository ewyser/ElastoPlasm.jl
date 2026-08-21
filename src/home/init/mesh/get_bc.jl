"""
    set_dirichlet(bc, xn, lower_upper_lim, dim; node::Symbol=:roller)

Flag Dirichlet boundary nodes in `bc` on dimension `dim`, at the coordinate `lower_upper_lim`.

# Arguments
- `bc`: Boolean matrix (ndim × nno) to update in-place.
- `xn`: Matrix of nodal coordinates.
- `lower_upper_lim`: Boundary coordinate value to match nodes against.
- `dim`: Dimension along which to apply the condition.
- `node::Symbol=:roller`: `:roller` flags only dimension `dim`; `:fixed` flags all dimensions.

# Returns
- `nothing`. Updates `bc` in-place.
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
    get_bc(xn::Matrix{T2}, solver::S) where {S<:AbstractSolver}

Compute the boundary condition matrix and boundary coordinates for a mesh, using boundary condition types specified in `solver.bcs.dirichlet`.

# Arguments
- `xn::Matrix{T2}`: Matrix of nodal coordinates (size: ndim × nno).
- `solver::S`: Solver instance (e.g. `ExplicitSolver`), must have `bcs.dirichlet` as a matrix of Symbols (e.g., `:roller`, `:fixed`).

# Returns
- `bc::Matrix{Bool}`: Boolean matrix (ndim × nno) indicating where boundary conditions are applied.
- `xB::Matrix{T2}`: Matrix (ndim × 2) of minimum and maximum boundary coordinates.

# Example
```julia
bc, xB = get_bc(xn, solver)
```
"""
function get_bc(xn::Matrix{T2},solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    xB = hcat(
        minimum(xn,dims=2),
        maximum(xn,dims=2)
    )
    ndim = size(xn,1)
    nno  = size(xn,2)
    bc   = zeros(Bool,ndim,nno)
    for dim ∈ 1:ndim
        for (k,limit) ∈ enumerate(xB[dim,:])
            set_dirichlet(bc,xn,limit,dim; node = solver.bcs.dirichlet[dim,k])
        end
    end
    return bc,xB
end