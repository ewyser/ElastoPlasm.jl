function set_roller_dirichlet(nno,xn,xB)
    bc = zeros(Bool,size(xB,1),nno[end])
    if size(xB,1) == 1

    elseif size(xB,1) == 2
        idx       = findall(x-> x ∈ xB[1,:],xn[1,:])
        idz       = findall(x-> x ∈ xB[2,:],xn[2,:])
        bc[1,idx].= true
        bc[2,idz].= true
    elseif size(xB,1) == 3
        idx = findall(x-> x ∈ xB[1,:],xn[1,:])
        idy = findall(x-> x ∈ xB[2,:],xn[2,:])
        idz = findall(x-> x ∈ xB[3,:],xn[3,:])
        bc[1,idx].= true
        bc[2,idy].= true
        bc[3,idz].= true
    end
    return bc
end
function set_roller_fixed(nno,xn,xB)
    bc = zeros(Bool,size(xB,1),nno[end])
    if size(xB,1) == 1

    elseif size(xB,1) == 2
        idx      = findall(x-> x ∈ xB[1,:],xn[1,:])
        idz      = findall(x-> x ∈ xB[2,:],xn[2,:])
        id       = unique(vcat(idx,idz))
        bc[1,id].= true
        bc[2,id].= true
    elseif size(xB,1) == 3
        idx = findall(x-> x ∈ xB[1,:],xn[1,:])
        idy = findall(x-> x ∈ xB[2,:],xn[2,:])
        idz = findall(x-> x ∈ xB[3,:],xn[3,:])
        id  = unique(vcat(idx,idy,idz))
        bc[1,id].= true
        bc[2,id].= true
        bc[3,id].= true
    end
    return bc
end

function get_bc(xn::Matrix{T2},nno::Vector{T1},ndim::T1; ghosts::Vector{T2}=[T2(0.0)], set::Symbol=:roller) where {T1,T2}
    l  = minimum(xn,dims=2).+ghosts
    L  = maximum(xn,dims=2).-ghosts

    xB = hcat(l,L)
    if set == :roller
        bc = set_roller_dirichlet(nno,xn,xB)
    elseif set == :fixed
        bc = set_fixed_dirichlet(nno,xn,xB)
    else
        throw(error("UnsupportedBC: $(set)"))
    end
    return bc,xB
end