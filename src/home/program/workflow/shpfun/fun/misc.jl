@kernel inbounds = true function Δnp_nd(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.prprt.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute delta functions
            for i ∈ 1:mesh.prprt.dim
                mpts.Δnp[nn,i,p] = mesh.x[i,no]-mpts.x[i,p]
            end
        end
    end
end
@views @kernel inbounds = true function Dij_nd(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.prprt.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute Dᵢⱼ tensor
            δ                = mesh.x[:,no].-mpts.x[:,p]
            mpts.Δnp[nn,:,p].= δ
            mpts.Dᵢⱼ[:,:,p] .= mpts.ϕ∂ϕ[nn,p,1].*(δ*δ')
        end
    end
end