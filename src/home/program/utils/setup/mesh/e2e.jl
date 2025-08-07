function e2e(ndim::T1,nel::Vector{T1},h::Vector{T2},instr) where {T1,T2}
    e2e  = spzeros(T1,nel[end],nel[end])
    nnel = ceil.(T1,instr[:nonloc][:ls]./h)
    gnum = collect(T1(1):nel[end])
    if ndim == 1
        iel  = 0
        for i ∈ 1:nel[1]#nelx
            iel = iel+1
            I   = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
            els = vec(gnum[I])         
            e2e[iel,els] = els
        end
    elseif ndim == 2
        gnum = reshape(gnum,Int(nel[2]),Int(nel[1]))
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
        gnum = reshape(gnum,Int(nel[3]),Int(nel[1]),Int(nel[2]))
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
