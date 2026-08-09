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
        J⁻¹ = T2(1.0)/mpts.J[p]
        mpts.s.σᵢ[p] = mpts.s.τᵢ[p] .* J⁻¹
    end   
end