"""
    flip_nd_n2p(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, dt::T2) where {T1,T2}

Update material point velocities and positions from mesh nodes (FLIP scheme, n-dim).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.

# Returns
- Updates material point fields in-place.
"""
@kernel inbounds = true function flip_nd_n2p(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        # flip update
        for dim ∈ 1:mesh.dim
            for nn ∈ 1:mesh.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                mpts.s.v[dim,p]+= dt*(mpts.ϕ∂ϕ[nn,p,1]*mesh.a[dim,no])
                mpts.x[dim,p]  += dt*(mpts.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no])
            end
            # find maximum velocity component over mps
            @atom mpts.vmax[dim] = max(mpts.vmax[dim],abs(mpts.s.v[dim,p]))
        end
    end  
end
"""
    pic_nd_n2p(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, dt::T2) where {T1,T2}

Update material point velocities and positions from mesh nodes (PIC scheme, n-dim).

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.

# Returns
- Updates material point fields in-place.
"""
@kernel inbounds = true function pic_nd_n2p(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mpts.nmp    
        # pic update
        for dim ∈ 1:mesh.dim
            δv = T2(0.0)
            for nn ∈ 1:mesh.nn
                no = mpts.p2n[nn,p]
                if iszero(no) continue end
                δv += mpts.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no]
            end
            mpts.s.v[dim,p] = δv 
            mpts.x[dim,p]  += dt*δv
            # find maximum velocity component over mps
            @atom mpts.vmax[dim] = max(mpts.vmax[dim],abs(mpts.s.v[dim,p]))
        end
    end  
end
"""
    n2p(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, dt::T2, instr::NamedTuple) where {T1,T2}

Map mesh node solution back to material points using the selected transfer kernel.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `dt::T2`: Time step.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function n2p(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # mapping to material point
    instr[:cairn][:mapsto][:map].n2p!(ndrange=mpts.nmp,mpts,mesh,dt);sync(CPU())
    return nothing
end