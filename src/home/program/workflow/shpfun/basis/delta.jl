@views @kernel inbounds = true function δ_1d(mp,mesh)
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mp.nmp
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            # compute delta functions
            mp.δnp[nn,1,p] = -(mp.x[p]-mesh.x[no])
        end
    end
end
@views @kernel inbounds = true function δ_2d(mp,mesh)
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mp.nmp
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            # compute delta functions
            mp.δnp[nn,1,p] = -(mp.x[1,p]-mesh.x[1,no])
            mp.δnp[nn,2,p] = -(mp.x[2,p]-mesh.x[2,no])
        end
    end
end
@views @kernel inbounds = true function δ_3d(mp,mesh)
    p = @index(Global)
    # calculate delta functions for tpic
    if p ≤ mp.nmp
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            # compute delta functions
            mp.δnp[nn,1,p] = -(mp.x[1,p]-mesh.x[1,no])
            mp.δnp[nn,2,p] = -(mp.x[2,p]-mesh.x[2,no])
            mp.δnp[nn,3,p] = -(mp.x[3,p]-mesh.x[3,no])
        end
    end
end