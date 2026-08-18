@kernel inbounds = true function shpfun!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D,NN}) where {T1,T2,D,NN,E,R}
    p = @index(Global)
    # cache shape values/gradients for all NN neighbors of particle p, once per timestep
    if p ≤ mpts.nmp
        Nbuf  = MVector{NN,T2}(undef)
        ∂Nbuf = MMatrix{NN,D,T2}(undef)
        for nn ∈ 1:NN
            _, N, ∂N     = eval_basis(mpts, mesh, basis, p, T1(nn))
            Nbuf[nn]     = N
            ∂Nbuf[nn,:] .= ∂N
        end
        basis.N[p]  = SVector{NN,T2}(Nbuf)
        basis.∂N[p] = SMatrix{NN,D,T2}(∂Nbuf)
    end
end

@views @kernel inbounds = true function Dij_nd(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D}) where {T1,T2,E,R,D}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N  = basis.N[p][nn]
            # compute Dᵢⱼ tensor
            δ                = mesh.x[no] .- mpts.x[p]
            mpts.Δnp[nn,:,p].= δ
            mpts.Dᵢⱼ[:,:,p] .= N.*(δ*δ')
        end
    end
end