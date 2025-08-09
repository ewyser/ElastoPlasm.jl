"""
    p2e2n_1d(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}) where {T1,T2}

Assign 1D material points to elements and nodes (topology kernel).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.

# Returns
- Updates connectivity fields in-place.
"""
@kernel inbounds = true function p2e2n_1d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # Compute element indices
        elx = clamp(fld(mpts.x[1,p]-mesh.x₀[1],mesh.h[1]),0,mesh.nel[1]-1)
        el  = round(T1,1+elx)
        #=
        elx = unsafe_trunc(T1,floor((mpts.x[1,p]-mesh.x₀[1])/mesh.h[1]))
        el  = round(T1,1+elx)  # +1 because Julia is 1-based 
        =#
        # Assign mpts-to-node connectivity
        for nn ∈ 1:mesh.nn
            mpts.p2n[nn,p] = mesh.e2n[nn,el]
        end
        # Store element index in mpts
        mpts.p2e[p] = el
    end
end
"""
    p2e2n_2d(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}) where {T1,T2}

Assign 2D material points to elements and nodes (topology kernel).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.

# Returns
- Updates connectivity fields in-place.
"""
@kernel inbounds = true function p2e2n_2d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # Compute element indices
        elx = clamp(fld(mpts.x[1,p]-mesh.x₀[1],mesh.h[1]),0,mesh.nel[1]-1)
        ely = clamp(fld(mpts.x[2,p]-mesh.x₀[2],mesh.h[2]),0,mesh.nel[2]-1)
        el  = round(T1,1+ely+mesh.nel[2]*elx)  # +1 because Julia is 1-based 
        #=
        elx = unsafe_trunc(T1,floor((mpts.x[1,p]-mesh.x₀[1])/mesh.h[1]))
        ely = unsafe_trunc(T1,floor((mpts.x[2,p]-mesh.x₀[2])/mesh.h[2]))
        el  = round(T1,1+ely+mesh.nel[2]*elx)  # +1 because Julia is 1-based 
        =#
        # Assign mpts-to-node connectivity
        for nn ∈ 1:mesh.nn
            mpts.p2n[nn,p] = mesh.e2n[nn,el]
        end
        # Store element index in mpts
        mpts.p2e[p] = el
    end
end
"""
    p2e2n_3d(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}) where {T1,T2}

Assign 3D material points to elements and nodes (topology kernel).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.

# Returns
- Updates connectivity fields in-place.
"""
@kernel inbounds = true function p2e2n_3d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp 
        # Compute element indices
        elx = clamp(fld(mpts.x[1,p]-mesh.x₀[1],mesh.h[1]),0,mesh.nel[1]-1)
        ely = clamp(fld(mpts.x[2,p]-mesh.x₀[2],mesh.h[2]),0,mesh.nel[2]-1)
        elz = clamp(fld(mpts.x[3,p]-mesh.x₀[3],mesh.h[3]),0,mesh.nel[3]-1)
        el  = round(T1,1+elz+mesh.nel[3]*elx+mesh.nel[3]*mesh.nel[1]*ely)
        #=
        elx = unsafe_trunc(T1,floor((mpts.x[1,p]-mesh.x₀[1])/mesh.h[1]))
        ely = unsafe_trunc(T1,floor((mpts.x[2,p]-mesh.x₀[2])/mesh.h[2]))
        elz = unsafe_trunc(T1,floor((mpts.x[3,p]-mesh.x₀[3])/mesh.h[3]))
        el  = round(T1,1+elz+mesh.nel[3]*elx+mesh.nel[3]*mesh.nel[1]*ely)
        =#
        # Assign mpts-to-node connectivity
        for nn ∈ 1:mesh.nn
            mpts.p2n[nn,p] = mesh.e2n[nn,el]
        end
        # Store element index in mpts
        mpts.p2e[p] = el
    end
end