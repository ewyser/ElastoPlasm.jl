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
        # σ = τ/J : scale the stored (p,dev) split of the Kirchhoff stress, giving a
        # CauchyStress whose `get_vector` equals the old `mpts.s.τᵢ[p] .* J⁻¹` Voigt product.
        mpts.s.σᵢ[p] = CauchyStress(mpts.s.τᵢ[p], J⁻¹)
    end   
end