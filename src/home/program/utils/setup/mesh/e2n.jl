"""
    e2n(ndim::T1, nno::Vector{T1}, nel::Vector{T1}, nn::T1) where {T1}

Construct the element-to-node connectivity matrix for a structured mesh.

# Arguments
- `ndim::T1`: Number of spatial dimensions.
- `nno::Vector{T1}`: Number of nodes in each direction and total.
- `nel::Vector{T1}`: Number of elements in each direction and total.
- `nn::T1`: Total number of nodes.

# Returns
- `e2n::Matrix{T1}`: Element-to-node connectivity matrix.

# Example
```julia
e2n_mat = e2n(2, [11, 11, 121], [10, 10, 100], 121)
```
"""
function e2n(ndim::T1,nno::Vector{T1},nel::Vector{T1},nn::T1) where {T1}
	iel,e2n =1,zeros(T1,nn,nel[end])
    gnum = collect(T1(1):nno[end])
    if ndim == 1
        for i0 ∈ 1:nel[1]#nelx
            nno = []
            for i ∈ -1:2
                try
                    push!(nno,gnum[i0+i])
                catch
                    push!(nno,0)
                end
            end
            e2n[:,iel].= nno
            iel        = iel+1
        end
    elseif ndim == 2
        gnum = reshape(gnum,Int(nno[2]),Int(nno[1]))
        for i0 ∈ 1:nel[1]#nelx
            for j0 ∈ 1:nel[2]#nelz
                nno = []
                for i ∈ -1:2
                    for j ∈ -1:2
                        try
                            push!(nno,gnum[j0+j,i0+i])
                        catch
                            push!(nno,0)
                        end
                    end
                end
                e2n[:,iel].= nno
                iel        = iel+1
            end
        end
    elseif ndim == 3
        gnum = reshape(gnum,Int(nno[3]),Int(nno[1]),Int(nno[2]))
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
                                    push!(nno,0)
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