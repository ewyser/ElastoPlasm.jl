@views @kernel inbounds = true function Dij_nd(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},basis::Basis{T1,D}) where {T1,T2,E,R,D}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = eval_basis(mpts, mesh, basis, p, nn)
            if iszero(no) continue end
            # compute Dᵢⱼ tensor
            δ                = mesh.x[no] .- mpts.x[p]
            mpts.Δnp[nn,:,p].= δ
            mpts.Dᵢⱼ[:,:,p] .= N.*(δ*δ')
        end
    end
end