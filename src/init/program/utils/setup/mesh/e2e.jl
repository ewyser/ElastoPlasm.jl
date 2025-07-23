function e2e(ndim,nno,nel,nn,h,instr)
    e2e  = spzeros(Int64,nel[end],nel[end])
    nnel = ceil.(Int,instr[:nonloc][:ls]./h)
    if ndim == 1
        gnum = collect(1:nel[end])
        iel  = 0
        for i ∈ 1:nel[1]#nelx
            iel = iel+1
            I   = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
            els = vec(gnum[I])         
            e2e[iel,els] = els
        end
    elseif ndim == 2
        gnum = reshape(1:nel[end],nel[2],nel[1])
        iel  = 0
        for i ∈ 1:nel[1]#nelx
            for j ∈ 1:nel[2]#nelz
                iel = iel+1
                I   = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
                J   = max(1,j-nnel[2]):min(nel[2],j+nnel[2])
                els = vec(gnum[J,I])         
                e2e[els,iel] = els
            end
        end
    elseif ndim == 3
        gnum = reshape(1:(nel[end]),nel[3],nel[1],nel[2])
        iel  = 0
        for k ∈ 1:nel[2] #nely
            for i ∈ 1:nel[1] #nelx
                for j ∈ 1:nel[3] #nelz
                    iel = iel+1
                    I   = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
                    J   = max(1,j-nnel[3]):min(nel[3],j+nnel[3])
                    K   = max(1,k-nnel[2]):min(nel[2],k+nnel[2])
                    els = vec(gnum[J,I,K])         
                    e2e[iel,els] = els
                end
            end
        end
    end
	return e2e
end
