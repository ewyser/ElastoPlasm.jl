@kernel inbounds = true function shpfun!(mpts::Point{T1,T2,D,CM},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D,NN}) where {T1,T2,D,NN,CM}
    p = @index(Global)
    # cache shape values/gradients for all NN neighbors of particle p, once per timestep
    if p ≤ mpts.nmp
        for nn ∈ 1:NN
            _, N, ∂N = eval_basis(mpts, mesh, basis, T1(p), T1(nn))
            basis.N[nn,p] = N
            for d ∈ 1:D
                basis.∂N[nn,d,p] = ∂N[d]
            end
        end
    end
end

@kernel inbounds = true function Dij_nd(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D}) where {T1,T2,D}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        Dᵢⱼ   = zero(eltype(mpts.Dᵢⱼ))
        nodes = basis.p2n[p]
        for nn ∈ 1:mesh.prprt.nn
            no  = nodes[nn]
            if iszero(no) continue end
            N   = basis.N[nn,p]
            # accumulate Dᵢⱼ tensor over all neighbors
            δ   = mesh.x[no] - mpts.x[p]
            Dᵢⱼ = Dᵢⱼ + N*(δ*δ')
        end
        mpts.Dᵢⱼ[p] = Dᵢⱼ
    end
end