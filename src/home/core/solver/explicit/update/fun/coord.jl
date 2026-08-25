# Generic: length(u) is a compile-time constant so the loop unrolls.
@kernel inbounds = true function coord(mpts::Point{T1,T2,D,E},mesh::Mesh{T1,T2,D}) where {T1,T2,D,E}
    p = @index(Global)
    if p ≤ mpts.nmp
        mpts.x[p] += mpts.s.u[p]
    end
end