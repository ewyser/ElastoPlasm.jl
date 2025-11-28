export get_glob, set_dispatch
"""
    get_glob(dims::Tuple,nrank::Int; zero::Bool=true)

Find optimal global matrix dimension(s) regarding the number of mpi processes given by `nrank`.

Returns a named tuple with:
- `dim`: tuple of dimension(s)
- `glob`: uninitialized Vector or Matrix of NamedTuples

Example:
---
```julia-repl
julia> set_glob_MPI(3)
(dim = (3,), glob = NamedTuple[#undef, #undef, #undef])

julia> 
```
"""
function get_glob(grid::Tuple,nrank::Int; zero::Bool=true)
    function find_optimal(nrank::Int)
        root = floor(Int, sqrt(nrank))
        for i ∈ root:-1:2
            if nrank % i == 0
                return (i, div(nrank, i))
            end
        end
        return (nrank,)
    end
    function build_glob(grid::Tuple,nprocs::Tuple; zero::Bool=true)
        ranks     = if zero collect(0:prod(nprocs)-1) else collect(1:prod(nprocs)) end
        Δnxy      = round.(Int,grid./nprocs)
        topology  = fill(MPI.PROC_NULL,nprocs.+2)
        dispatch  = Array{Dict,length(nprocs)}(undef,nprocs)
        neighbors = Array{Dict,length(nprocs)}(undef,nprocs)
        indices   = Array{Dict,length(nprocs)}(undef,nprocs)
        if length(nprocs) == 1
            topology[2:end-1].= ranks
            idx,idy   = collect(1:Δnxy[1]),collect(1:grid[2])
            for i ∈ 2:size(topology,1)-1
                dispatch[i-1]  = Dict(
                    :me => topology[i],
                )
                neighbors[i-1] = Dict(
                    :dir => Dict(
                        :x=> Dict(
                            :me => topology[i],
                            :send  => topology[i-1],
                            :recv => topology[i+1],
                        ),
                    ),
                )
                indices[i-1] = Dict(
                    :me => topology[i],
                    :i => (idx[1]:idx[end]),
                    :j => (idy[1]:idy[end]),
                )
                idx = unique(min.(idx.+Δnxy[1],grid[1]))
            end
        elseif length(nprocs) == 2
            topology[2:end-1,2:end-1].= reshape(ranks,nprocs[2],nprocs[1])'
            idx = collect(1:Δnxy[1])
            for i ∈ 2:size(topology,1)-1
                idy = collect(1:Δnxy[2])
                for j ∈ 2:size(topology,2)-1
                    dispatch[i-1,j-1]  = Dict(
                        :me => topology[i,j],
                    )
                    neighbors[i-1,j-1] = Dict(
                        :dir => Dict(
                            :x=> Dict(
                                :me => topology[i,j],
                                :send  => topology[i-1,j],
                                :recv => topology[i+1,j],
                            ),
                            :y=> Dict(
                                :me => topology[i,j],
                                :send  => topology[i,j-1],
                                :recv => topology[i,j+1],
                            ),
                        ),
                    )
                    indices[i-1,j-1] = Dict(
                        :me => topology[i,j],
                        :glob => Dict(
                            :i => (idx[1]:idx[end]),
                            :j => (idy[1]:idy[end]),
                        ),
                        :loc => Dict(
                            :i => (1:length(idx)),
                            :j => (1:length(idy)),
                        ),
                        :i => (idx[1]:idx[end]),
                        :j => (idy[1]:idy[end]),
                    )
                    idy = unique(min.(idy.+Δnxy[2],grid[2]))
                end
                idx = unique(min.(idx.+Δnxy[1],grid[1]))
            end
        end
        return Dict(
            :tplgy   => if length(nprocs) == 1 topology[2:end-1] else topology[2:end-1, 2:end-1] end,
            :neighb  => vec(permutedims(neighbors)),
            :index   => vec(permutedims(indices)),
            :globdim => Int.((grid[1],grid[2])),
        )
    end
    nprocs = find_optimal(nrank)
    # initialize 
    return build_glob(grid,nprocs; zero=zero)
end 
"""
    set_dispatch(startup::String)

Find optimal global matrix dimension(s) regarding the number of mpi processes given by `nrank`.

Returns a named tuple with:
- `dim`: tuple of dimension(s)
- `glob`: uninitialized Vector or Matrix of NamedTuples

Example:
---
```julia-repl
julia> set_glob_MPI(3)
(dim = (3,), glob = NamedTuple[#undef, #undef, #undef])

julia> 
```
"""
function set_dispatch(startup::String)
    mpi = load(startup)["config"].mpi
    if mpi[1].status
        nproc  = length(mpi)
        try CUDA.functional()
            expr = "using cORIUm, CUDA;"
        catch 
            expr = "using cORIUm;"
        end
        if Sys.iswindows()
            expr = "$expr cORES!(raw\"$startup\")"
        else
            expr = "$expr cORES!(\"$startup\")"
        end
        mpiexec(cmd -> run(`$cmd -n $(nproc) julia --project=. -e $expr `))
        return load(startup)["ic"]::String,load(startup)["h5file"]::String
    else
        id     = floor(Int,rand()*length(mpi))+1
        mpi    = mpi[id]
        ic,dat = cORE!(mpi,startup)
        return ic::String, dat::String
    end
end