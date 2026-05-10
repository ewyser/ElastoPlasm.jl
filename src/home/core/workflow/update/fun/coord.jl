@kernel inbounds = true function coord(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D}) where {T1,T2,D<:TwoDimension,E,R}
    p = @index(Global)
    if p ≤ mpts.nmp 
        mpts.x[1,p]+= mpts.s.u[1,p]
        mpts.x[2,p]+= mpts.s.u[2,p]
    end
end
@kernel inbounds = true function coord(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D}) where {T1,T2,D<:ThreeDimension,E,R}
    p = @index(Global)
    if p ≤ mpts.nmp 
        mpts.x[1,p]+= mpts.s.u[1,p]
        mpts.x[2,p]+= mpts.s.u[2,p]
        mpts.x[3,p]+= mpts.s.u[3,p]
    end
end