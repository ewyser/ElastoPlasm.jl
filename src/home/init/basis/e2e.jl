"""
    e2e(nel::Vector{T1}, h::Vector{T2}, solver::S) where {T1,T2,S<:AbstractSolver}

Construct the element-to-element connectivity matrix for a structured mesh, in any dimension:
column `el` holds, as row indices, every element within `ceil(nonloc.ls/h)` elements of `el` along
each axis (itself included). Only the sparsity pattern is read (`nzrange`/`rowvals` in
`nonlocal.jl`). Built in one `sparse(I, J, V)` call: filling a `SparseMatrixCSC` entry by entry
costs `O(nnz)` per insertion.

# Arguments
- `nel::Vector{T1}`: Number of elements in each direction, then the total.
- `h::Vector{T2}`: Element size in each direction.
- `solver::S`: Solver instance (e.g. `ExplicitSolver`), must have `nonloc.ls` (nonlocal support length scale).

# Returns
- `e2e::SparseMatrixCSC{T1,T1}`: Element-to-element connectivity matrix.

# Example
```julia
e2e_mat = e2e([10, 10, 100], [0.1, 0.1], solver)
```
"""
function e2e(nel::Vector{T1},h::Vector{T2},solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    n    = ntuple(d -> Int(nel[d]), D)
    r    = CartesianIndex(ntuple(d -> ceil(Int, solver.nonloc.ls / h[d]), D))
    C, L = CartesianIndices(n), LinearIndices(n)
    I, J = sizehint!(T1[], prod(2 .* Tuple(r) .+ 1) * length(C)), sizehint!(T1[], prod(2 .* Tuple(r) .+ 1) * length(C))
    for c ∈ C, c′ ∈ max(c - r, first(C)):min(c + r, last(C))
        push!(I, T1(L[c′])); push!(J, T1(L[c]))
    end
    return sparse(I, J, I, T1(nel[end]), T1(nel[end]))
end
