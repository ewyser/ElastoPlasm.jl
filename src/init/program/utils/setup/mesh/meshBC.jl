function meshBCs(xn,h,nno,nD;drift=0.0)
    l = minimum(xn,dims=1).+drift
    L = maximum(xn,dims=1).-drift
    if nD == 1
        xB  = [minimum(xn)+drift*h,maximum(xn)-drift*h]
        bcX = ones(Float64,nno[end])
        bcX[1]  = 0.0
        bcX[end]= 0.0
        bc   = bcX
    elseif nD == 2
        xB  = vcat([l[1],L[1]],[l[2],Inf])
        bcx = findall(x-> x ∈ xB[1:2],xn[:,1])
        bcz = findall(x-> x ∈ xB[3:4],xn[:,2])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        bc   = hcat(bcX,bcZ)
    elseif nD == 3
        xB  = vcat([l[1],L[1]],[l[2],L[2]],[l[3],Inf])
        bcx = findall(x-> x ∈ xB[1:2],xn[:,1])
        bcy = findall(x-> x ∈ xB[3:4],xn[:,2])
        bcz = findall(x-> x ∈ xB[5:6],xn[:,3])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcY = ones(nno[nD+1],1)
        bcY[bcy] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        bc   = hcat(bcX,bcY,bcZ)
    end
    return bc,xB
end