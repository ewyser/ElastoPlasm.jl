function e2n(ndim,nno,nel,nn)
	iel,e2n =1,zeros(Int64,nn,nel[end])
    if ndim == 1
        gnum = collect(1:nno[end])
        for i0 ∈ 1:nel[1]#nelx
            nno = []
            for i ∈ -1:2
                try
                    push!(nno,gnum[i0+i])
                catch
                    push!(nno,-404)
                end
            end
            e2n[:,iel].= nno
            iel        = iel+1
        end
    elseif ndim == 2
        gnum = reshape(1:(nno[end]),nno[2],nno[1])
        for i0 ∈ 1:nel[1]#nelx
            for j0 ∈ 1:nel[2]#nelz
                nno = []
                for i ∈ -1:2
                    for j ∈ -1:2
                        try
                            push!(nno,gnum[j0+j,i0+i])
                        catch
                            push!(nno,-404)
                        end
                    end
                end
                e2n[:,iel].= nno
                iel        = iel+1
            end
        end
    elseif ndim == 3
        gnum = reshape(1:(nno[end]),nno[3],nno[1],nno[2])
        for k0 ∈ 1:nel[2]#nely
            for i0 ∈ 1:nel[1]#nelx
                for j0 ∈ 1:nel[3]#nelz gnum[j0-1,i0-1,k0-1]
                    nno = []
                    for k ∈ -1:2
                        for i ∈ -1:2
                            for j ∈ -1:2
                                try
                                    push!(nno,gnum[j0+j,i0+i,k0+k])
                                catch
                                    push!(nno,-404)
                                end
                            end
                        end
                    end
                    e2n[:,iel].= nno
                    iel        = iel+1
                end
            end
        end
    end
	return e2n
end