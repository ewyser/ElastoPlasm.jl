@kernel inbounds = true function mass_assembly(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        for nn ∈ 1:stencil_size(basis)
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            @atom mesh.s.m[no] += basis.N[nn,p] * mpts.s.m[p]
        end
    end
end
