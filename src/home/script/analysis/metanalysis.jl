using Suppressor

export metanalysis

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
function metanalysis(L,nel; fid::String=first(splitext(basename(@__FILE__))),)
    # Setting up cases for metanalysis
    deforms,trfrs = ["finite", "infinitesimal"],["std", "tpic", "apic"]
    sim,nsim      = 1,length(deforms)*length(trfrs)    
    
    # Preparing mesh and material points for metanalysis
    @info "Running metanalysis for $(length(L))d slump problem"
    prog = make_progress(nsim; desc="Preparing simulation(s)")
    for (i,deform) ∈ enumerate(deforms)
        for (j,trsfr) ∈ enumerate(trfrs)
            fwrk  = (;
                deform = deform,
                trsfr = trsfr,
                C_pf = 1.0, 
                musl = true,
                locking = true,
                damping = 0.1
            )
            plot  = (;
                status = false,
                freq   = 1.0,
                dpi    = 500,
                what   = [(;mpts=(name="epII",cblim=(0.0,1.5)),),],
            )

            @suppress begin
                _ = ic_slump(L,nel;fid="$fid/sim_$sim",fwrk=fwrk,plot=plot);
            end            
            next!(prog;desc="Preparing simulation(s) $sim/$nsim...")
            sim += 1
        end
    end
    
    # Get all subdirectories in the meta folder
    meta_path = joinpath(info.sys.out, fid)
    sims = filter(d -> isdir(joinpath(meta_path, d)), readdir(meta_path))
    @info "Available simulations for metanalysis in $(trunc_path(meta_path)):"
    outputs = []
    prog = make_progress(nsim; desc="Executing simulation(s)")
    for (k, sim) ∈ enumerate(sims)
        sim_path = joinpath(meta_path, sim)
        for file in filter(f -> endswith(f, ".jld2"), readdir(sim_path))
            @suppress begin
                out  = elastoplasm!(joinpath(sim_path,file); workflow = [elastodynamic!,elastoplastic!]);
                push!(outputs, out)
            end  
        end
        next!(prog;desc="Executing simulation(s) $k/$nsim...")
    end
    println(outputs)
    #TODO: add metanalysis post-processing here
    return nothing
end