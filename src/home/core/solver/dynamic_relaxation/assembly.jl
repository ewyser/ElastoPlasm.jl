@kernel inbounds = true function mass_assembly(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},g::Vector{T2}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        # buffering 
        ms = mpts.s.ρ[p]*mpts.Ω[p]
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            # accumulation
            @atom mesh.s.m[no] += N * ms
        end
    end
end