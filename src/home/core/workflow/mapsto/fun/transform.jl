"""
    transform(mpts::Point{T1,T2}) where {T1,T2}

Transform Kirchhoff to Cauchy stress at material points.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.

# Returns
- Updates stress in-place.
"""
@views @kernel inbounds = true function transform(mpts::Point{T1,T2}) where {T1,T2}
    p = @index(Global)
    # deformation framework dispatcher
    if p ≤ mpts.nmp 
        J = mpts.J[p]
        for c ∈ 1:size(mpts.s.σᵢ,1)
            mpts.s.σᵢ[c,p] = mpts.s.τᵢ[c,p]/J
        end
    end   
end