"""
    e2e(nel, h, solver) -> Adjacency

The elements near each element: column `el` holds, as row indices, every element within
`ceil(nonloc.ls/h)` elements of `el` along each axis (itself included), for the nonlocal
regularisation (`nonlocal.jl` reads only the pattern, `nzrange`/`rowvals`). Built in one `sparse(I, J, V)` call,
since inserting entries one by one costs `O(nnz)` each, and kept as its pattern, an `Adjacency` (which a device holds).
"""
function e2e(nel::Vector{T1},h::Vector{T2},solver::S) where {T1,T2,D,S<:Solver{T1,T2,D}}
    n    = ntuple(d -> Int(nel[d]), D)
    r    = CartesianIndex(ntuple(d -> ceil(Int, solver.nonloc.ls / h[d]), D))
    C, L = CartesianIndices(n), LinearIndices(n)
    I, J = sizehint!(T1[], prod(2 .* Tuple(r) .+ 1) * length(C)), sizehint!(T1[], prod(2 .* Tuple(r) .+ 1) * length(C))
    for c ∈ C, c′ ∈ max(c - r, first(C)):min(c + r, last(C))
        push!(I, T1(L[c′])); push!(J, T1(L[c]))
    end
    A = sparse(I, J, I, T1(nel[end]), T1(nel[end]))
    return Adjacency(A.colptr, A.rowval)
end
