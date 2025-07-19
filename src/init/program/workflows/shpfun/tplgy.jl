@kernel inbounds = true function p2e2n1D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.p2e[p] = cld(mpD.x[p,1]-meD.x₀[1],meD.h[1])
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
@kernel inbounds = true function p2e2n2D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        # Compute element indices
        elx = clamp(fld(mpD.x[p,1]-meD.x₀[1],meD.h[1]),0,meD.nel[1]-1)
        ely = clamp(fld(mpD.x[p,2]-meD.x₀[2],meD.h[2]),0,meD.nel[2]-1)

        # Compute flattened element index (column-major ordering)
        mpD.p2e[p] = 1+ely+meD.nel[2]*elx  # +1 because Julia is 1-based

        # Assign connectivity
        for nn ∈ 1:meD.nn
            mpD.p2n[nn, p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
@kernel inbounds = true function p2e2n3D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        # Compute element indices
        elx = clamp(fld(mpD.x[p,1]-meD.x₀[1],meD.h[1]), 0, meD.nel[1]-1)
        ely = clamp(fld(mpD.x[p,2]-meD.x₀[2],meD.h[2]), 0, meD.nel[2]-1)
        elz = clamp(fld(mpD.x[p,3]-meD.x₀[3],meD.h[3]), 0, meD.nel[3]-1)

        # Compute flattened element index (column-major ordering)
        mpD.p2e[p] = 1+elz+meD.nel[3]*elx+meD.nel[3]*meD.nel[1]*ely

        # Assign connectivity
        for nn in 1:meD.nn
            mpD.p2n[nn, p] = meD.e2n[nn, mpD.p2e[p]]
        end
    end
end