@kernel inbounds = true function p2e2n_1d(mp::Point{T1,T2},mesh) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp 
        # Compute element indices
        el = round(Int,cld(mp.x[p,1]-mesh.x₀[1],mesh.h[1]))
        # Assign mp-to-node connectivity
        for nn ∈ 1:mesh.nn
            mp.p2n[nn,p] = mesh.e2n[nn,el]
        end
        # Store element index in mp
        mp.p2e[p] = el
    end
end
@kernel inbounds = true function p2e2n_2d(mp::Point{T1,T2},mesh) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp 
        # Compute element indices
        elx = clamp(fld(mp.x[1,p]-mesh.x₀[1],mesh.h[1]),0,mesh.nel[1]-1)
        ely = clamp(fld(mp.x[2,p]-mesh.x₀[2],mesh.h[2]),0,mesh.nel[2]-1)
        el  = round(T1,1+ely+mesh.nel[2]*elx)  # +1 because Julia is 1-based 
        # Assign mp-to-node connectivity
        for nn ∈ 1:mesh.nn
            mp.p2n[nn,p] = mesh.e2n[nn,el]
        end
        # Store element index in mp
        mp.p2e[p] = el
    end
end
@kernel inbounds = true function p2e2n_3d(mp::Point{T1,T2},mesh) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp 
        # Compute element indices
        elx = clamp(fld(mp.x[1,p]-mesh.x₀[1],mesh.h[1]),0,mesh.nel[1]-1)
        ely = clamp(fld(mp.x[2,p]-mesh.x₀[2],mesh.h[2]),0,mesh.nel[2]-1)
        elz = clamp(fld(mp.x[3,p]-mesh.x₀[3],mesh.h[3]),0,mesh.nel[3]-1)
        el  = round(T1,1+elz+mesh.nel[3]*elx+mesh.nel[3]*mesh.nel[1]*ely)
        # Assign mp-to-node connectivity
        for nn ∈ 1:mesh.nn
            mp.p2n[nn,p] = mesh.e2n[nn,el]
        end
        # Store element index in mp
        mp.p2e[p] = el
    end
end