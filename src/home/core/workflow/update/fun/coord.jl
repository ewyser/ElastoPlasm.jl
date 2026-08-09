@kernel inbounds = true function coord(mpts::Point{T1,T2,D,<:AbstractBasis,E,R},mesh::Mesh{T1,T2,D}) where {T1,T2,D<:TwoDimension,E,R}
    p = @index(Global)
    if p ≤ mpts.nmp
        u = mpts.s.u[p]
        mpts.x[1,p] += u[1]
        mpts.x[2,p] += u[2]
    end
end
@kernel inbounds = true function coord(mpts::Point{T1,T2,D,<:AbstractBasis,E,R},mesh::Mesh{T1,T2,D}) where {T1,T2,D<:ThreeDimension,E,R}
    p = @index(Global)
    if p ≤ mpts.nmp
        u = mpts.s.u[p]
        mpts.x[1,p] += u[1]
        mpts.x[2,p] += u[2]
        mpts.x[3,p] += u[3]
    end
end