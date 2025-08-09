"""
    set_roller_dirichlet(nno, xn, xB)

Set roller Dirichlet boundary conditions for a mesh.

# Arguments
- `nno`: Number of nodes in each direction and total.
- `xn`: Matrix of nodal coordinates.
- `xB`: Boundary coordinates.

# Returns
- `bc::Matrix{Bool}`: Boolean matrix indicating roller boundary conditions.
"""
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
"""
    set_roller_fixed(nno, xn, xB)

Set roller fixed boundary conditions for a mesh.

# Arguments
- `nno`: Number of nodes in each direction and total.
- `xn`: Matrix of nodal coordinates.
- `xB`: Boundary coordinates.

# Returns
- `bc::Matrix{Bool}`: Boolean matrix indicating fixed boundary conditions.
"""
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

"""
    get_bc(xn::Matrix{T2}, nno::Vector{T1}, ndim::T1; ghosts::Vector{T2}=[T2(0.0)], set::Symbol=:roller) where {T1,T2}

Compute boundary condition matrix and boundary coordinates for a mesh.

# Arguments
- `xn::Matrix{T2}`: Matrix of nodal coordinates.
- `nno::Vector{T1}`: Number of nodes in each direction and total.
- `ndim::T1`: Number of spatial dimensions.
- `ghosts::Vector{T2}`: Optional, size of ghost cells to add at boundaries (default: `[T2(0.0)]`).
- `set::Symbol`: Type of boundary condition to set (`:roller` or `:fixed`).

# Returns
- `bc::Matrix{Bool}`: Boolean matrix indicating boundary conditions.
- `xB::Matrix{T2}`: Boundary coordinates.
"""
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