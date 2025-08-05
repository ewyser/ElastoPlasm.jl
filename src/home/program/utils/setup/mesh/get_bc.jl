function get_bc(xn,h,nno,ndim;ghosts=0.0)
    l = minimum(xn,dims=2).+ghosts
    L = maximum(xn,dims=2).-ghosts
    bc = ones(ndim,nno[end])
    if ndim == 1
        xB  = vcat(l,L)
        bcX = ones(Float64,nno[end])
        bcX[1]  = 0.0
        bcX[end]= 0.0
        bc   = bcX
    elseif ndim == 2
        xB  = vcat([l[1],L[1]],[l[2],L[2]])
        idx = findall(x-> x ∈ xB[1:2],xn[1,:])
        idz = findall(x-> x ∈ xB[3:4],xn[2,:])
        for k ∈ 1:nno[end]
            if k ∈ idx
                bc[1,k] = 0.0
            end
            if k ∈ idz
                bc[2,k] = 0.0
            end
        end
    elseif ndim == 3
        xB  = vcat([l[1],L[1]],[l[2],L[2]],[l[3],L[3]])
        idx = findall(x-> x ∈ xB[1:2],xn[1,:])
        idy = findall(x-> x ∈ xB[3:4],xn[2,:])
        idz = findall(x-> x ∈ xB[5:6],xn[3,:])
        for k ∈ 1:nno[end]
            if k ∈ idx
                bc[1,k] = 0.0
            end
            if k ∈ idy
                bc[2,k] = 0.0
            end
            if k ∈ idz
                bc[3,k] = 0.0
            end
        end
    end
    return bc,xB
end