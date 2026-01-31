export get_glob, set_dispatch

# Helper Functions
# ==============================================================================

"""
Find optimal 2D decomposition for nrank processes (returns (nrank,) if prime).
"""
function find_optimal_decomposition(nrank::Int)
    root = floor(Int, sqrt(nrank))
    for i ∈ root:-1:2
        if nrank % i == 0
            return (i, div(nrank, i))
        end
    end
    return (nrank,)
end

"""
Compute next index range, advancing by delta and clamping to grid limit.
"""
function next_indices(current::Vector, delta::Int, limit::Int)
    return unique(min.(current .+ delta, limit))
end

"""
Create neighbor info for a process at position (i,j) in topology.
"""
function make_neighbors(topology, i::Int, j::Int=0)
    if j == 0  # 1D case
        return Dict(
            :dir => Dict(
                :x => Dict(:me => topology[i], :send => topology[i-1], :recv => topology[i+1])
            )
        )
    else  # 2D case
        return Dict(
            :dir => Dict(
                :x => Dict(:me => topology[i,j], :send => topology[i-1,j], :recv => topology[i+1,j]),
                :y => Dict(:me => topology[i,j], :send => topology[i,j-1], :recv => topology[i,j+1])
            )
        )
    end
end

"""
Create index ranges for a process.
"""
function make_indices(rank::Int, x_range::UnitRange, y_range::UnitRange, is_2d::Bool=false)
    if !is_2d  # 1D case
        return Dict(
            :me => rank,
            :glob => Dict(:i => x_range, :j => y_range)
        )
    else  # 2D case
        return Dict(
            :me => rank,
            :glob => Dict(:i => x_range, :j => y_range),
            :loc => Dict(:i => 1:length(x_range), :j => 1:length(y_range))
        )
    end
end

# Main Functions
# ==============================================================================

"""
    get_glob(grid::Tuple, nrank::Int; zero::Bool=true) -> Dict

Decompose grid across nrank processes with optimal 2D distribution.

# Arguments
- `grid::Tuple`: Grid dimensions (nx, ny)
- `nrank::Int`: Number of processes
- `zero::Bool=true`: Use 0-based rank indexing

# Returns
- Dict with :tplgy (topology), :neighb (neighbors), :index (ranges), :globdim

# Example
```julia
decomp = get_glob((100, 50), 4)  # Splits into 2×2 process grid
```
"""
function get_glob(grid::Tuple, nrank::Int; zero::Bool=true)
    nprocs = find_optimal_decomposition(nrank)
    
    # Setup: rank numbering and grid partitioning
    ranks = zero ? (0:prod(nprocs)-1) : (1:prod(nprocs))
    cells_per_proc = round.(Int, grid ./ nprocs)
    PROC_NULL = isdefined(Main, :MPI) ? Main.MPI.PROC_NULL : -1
    
    # Build process topology (padded for MPI boundary communication)
    topology = fill(PROC_NULL, nprocs .+ 2)
    if length(nprocs) == 1
        topology[2:end-1] .= collect(ranks)
    else
        topology[2:end-1, 2:end-1] .= reshape(collect(ranks), nprocs[2], nprocs[1])'
    end
    
    # Allocate result storage
    n_neighbors = prod(nprocs)
    neighbors = Vector{Dict}(undef, n_neighbors)
    indices = Vector{Dict}(undef, n_neighbors)
    
    # Build decomposition
    if length(nprocs) == 1
        # 1D: Split along x-axis
        x_start = 1
        for (proc_id, i) in enumerate(2:size(topology,1)-1)
            x_end = min(x_start + cells_per_proc[1] - 1, grid[1])
            neighbors[proc_id] = make_neighbors(topology, i)
            indices[proc_id] = make_indices(topology[i], x_start:x_end, 1:grid[2], false)
            x_start = x_end + 1
        end
    else
        # 2D: Split along both axes
        x_start = 1
        proc_id = 0
        for i ∈ 2:size(topology,1)-1
            x_end = min(x_start + cells_per_proc[1] - 1, grid[1])
            y_start = 1
            for j ∈ 2:size(topology,2)-1
                y_end = min(y_start + cells_per_proc[2] - 1, grid[2])
                proc_id += 1
                neighbors[proc_id] = make_neighbors(topology, i, j)
                indices[proc_id] = make_indices(topology[i,j], x_start:x_end, y_start:y_end, true)
                y_start = y_end + 1
            end
            x_start = x_end + 1
        end
    end
    
    return Dict(
        :tplgy => length(nprocs) == 1 ? topology[2:end-1] : topology[2:end-1, 2:end-1],
        :neighb => neighbors,
        :index => indices,
        :globdim => Int.((grid[1], grid[2]))
    )
end