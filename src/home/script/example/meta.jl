export meta

"""
    ic_slump(L::Vector{Float64}, nel::Vector{Int64}; fid::String=..., kwargs...) -> NamedTuple, NamedTuple

Initializes the mesh, material points, and simulation configuration for a slump test.

# Arguments
- `L::Vector{Float64}`: Domain dimensions.
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration.

# Returns
- `(ic, cfg)`: Two named tuples containing mesh/material point/compression info (`ic`) and instructions/paths (`cfg`).

# Example
```julia
ic, cfg = ic_slump([64.0, 16.0], [40, 10]; fid="run1")
```
"""
function meta(L,nel; fid::String=first(splitext(basename(@__FILE__))),)
    @info "Running metanalysis for $(length(L))d slump problem"
    for (sim,deform) ∈ enumerate(["finite", "infinitesimal"])
        fwrk  = (;
            deform = deform,
            trsfr = "std",
            C_pf = 1.0, 
            musl = true,
            locking = true,
            damping = 0.1
        )
        jld2 = ic_slump(L,nel;fid="$fid/sim_$sim",fwrk=fwrk);
    end

    # Get all subdirectories in the meta folder
    meta_path = joinpath(info.sys.out, fid)
    sims = filter(d -> isdir(joinpath(meta_path, d)), readdir(meta_path))
    @info "Available simulations for metanalysis in $meta_path:"

    outputs = []
    for sim ∈ sims
        sim_path = joinpath(meta_path, sim)
        for file in filter(f -> endswith(f, ".jld2"), readdir(sim_path))
            jld2 = joinpath(sim_path, file)
            @info " Found simulation file: $file"
            out  = elastoplasm!(jld2; workflow = [elastodynamic!,elastoplastic!]);
            push!(outputs, out)
        end
    end
    println(outputs)
    #TODO: add metanalysis post-processing here
    return nothing
end