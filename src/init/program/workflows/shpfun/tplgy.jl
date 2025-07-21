@kernel inbounds = true function p2e2n_1d(mp,mesh)
    p = @index(Global)
    if p≤mp.nmp 
        mp.p2e[p] = cld(mp.x[p,1]-mesh.x₀[1],mesh.h[1])
        for nn ∈ 1:mesh.nn
            mp.p2n[nn,p] = mesh.e2n[nn,mp.p2e[p]]
        end
    end
end
@kernel inbounds = true function p2e2n_2d(mp,mesh)
    p = @index(Global)
    if p≤mp.nmp 
        # Compute element indices
        elx = clamp(fld(mp.x[1,p]-mesh.x₀[1],mesh.h[1]),0,mesh.nel[1]-1)
        ely = clamp(fld(mp.x[2,p]-mesh.x₀[2],mesh.h[2]),0,mesh.nel[2]-1)

        # Compute flattened element index (column-major ordering)
        mp.p2e[p] = 1+ely+mesh.nel[2]*elx  # +1 because Julia is 1-based

        # Assign connectivity
        for nn ∈ 1:mesh.nn
            mp.p2n[nn,p] = mesh.e2n[nn,mp.p2e[p]]
        end
    end
end
@kernel inbounds = true function p2e2n_3d(mp,mesh)
    p = @index(Global)
    if p≤mp.nmp 
        # Compute element indices
        elx = clamp(fld(mp.x[1,p]-mesh.x₀[1],mesh.h[1]),0,mesh.nel[1]-1)
        ely = clamp(fld(mp.x[2,p]-mesh.x₀[2],mesh.h[2]),0,mesh.nel[2]-1)
        elz = clamp(fld(mp.x[3,p]-mesh.x₀[3],mesh.h[3]),0,mesh.nel[3]-1)

        # Compute flattened element index (column-major ordering)
        mp.p2e[p] = 1+elz+mesh.nel[3]*elx+mesh.nel[3]*mesh.nel[1]*ely

        # Assign connectivity
        for nn ∈ 1:mesh.nn
            mp.p2n[nn,p] = mesh.e2n[nn,mp.p2e[p]]
        end
    end
end