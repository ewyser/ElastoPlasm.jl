"""
    set_dirichlet(bc, xn, limit, dim; node=:roller)

Flag the nodes of `xn` lying at `limit` along axis `dim` in `bc` (axes × nodes): `:roller` blocks
the motion along `dim` only, `:fixed` along every axis, `:free` nothing.
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
    get_bc(xn, solver) -> (bc, xB)

The Dirichlet flags of the nodes `xn` (axes × nodes, `true` where the motion along that axis is
blocked) from `solver.bcs.dirichlet` (one `[lower upper]` row per axis), and the domain bounds `xB`
(axes × 2) the nodes span.
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