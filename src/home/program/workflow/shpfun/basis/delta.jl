@views @kernel inbounds = true function δ_1d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute delta functions
            mpts.δnp[nn,1,p] = -(mpts.x[p]-mesh.x[no])
        end
    end
end
@views @kernel inbounds = true function δ_2d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute delta functions
            mpts.δnp[nn,1,p] = -(mpts.x[1,p]-mesh.x[1,no])
            mpts.δnp[nn,2,p] = -(mpts.x[2,p]-mesh.x[2,no])
        end
    end
end
@views @kernel inbounds = true function δ_3d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute delta functions
            mpts.δnp[nn,1,p] = -(mpts.x[1,p]-mesh.x[1,no])
            mpts.δnp[nn,2,p] = -(mpts.x[2,p]-mesh.x[2,no])
            mpts.δnp[nn,3,p] = -(mpts.x[3,p]-mesh.x[3,no])
        end
    end
end