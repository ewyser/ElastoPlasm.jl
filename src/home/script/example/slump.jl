export slump,ic_slump 

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
function ic_slump(L::Vector{Float64},nel::Vector{Int64}; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d slump problem"
    # init & kwargs
    instr   = kwargser(:instr,kwargs;dim=length(L))
    paths   = set_paths(fid,info.sys.out;interactive=false)      
    # mesh & mp initial conditions
    ni      = 2  
    mesh    = setup_mesh(nel,L,instr)    
    cmpr    = setup_cmpr(mesh.dim,instr)                       
    mp      = setup_mps(mesh,cmpr;define=inislump(mesh,cmpr,ni,instr))
    # time parameters
    time = (; T = 15.0, te = 10.0, tg = 15.0/1.5,)  
    return (;mesh,mp,cmpr,time),(;instr,paths)
end

"""
    slump(ic::NamedTuple, cfg::NamedTuple) -> Bool

Runs the explicit solution workflow for the slump problem, including simulation and postprocessing.

# Arguments
- `ic::NamedTuple`: Initial mesh/material/compression configuration.
- `cfg::NamedTuple`: Simulation instructions and output paths.

# Returns
- `Bool`: Returns `true` if the simulation and postprocessing complete successfully.

# Example
```julia
success = slump(ic, cfg)
```
"""
function slump(ic::NamedTuple,cfg::NamedTuple;)
    @info "Explicit solution to slump problem";config_plot()
    # extract mesh, mp, cmpr, instr, paths
    mesh,mp,cmpr = deepcopy(ic[:mesh]  ), deepcopy(ic[:mp]    ), deepcopy(ic[:cmpr])
    time         = deepcopy(ic[:time]  )
    instr,paths  = deepcopy(cfg[:instr]), deepcopy(cfg[:paths])
    # plot initial cohesion field
    plotcoh(mp,cmpr,paths)                                                  
    # action
    out  = plasming!(mp,mesh,cmpr,time,instr)
    # postprocessing
    @info "Fig(s) saved at $(paths[:plot])"
    path =joinpath(paths[:plot],"$(mesh.dim)d_$(mp.nmp)_$(mesh.nel[end])_$(join(instr[:plot][:what]))_$(instr[:basis][:which])_$(instr[:fwrk][:deform])_$(instr[:fwrk][:trsfr])_$(instr[:fwrk][:locking])_$(cmpr[:cmType])_$(instr[:perf])_$(first(instr[:nonloc])).png")
    savefig(path)
    msg("(✓) Done! exiting...\n")
    return sucess=true
end

#=
    L,nel  = [64.1584,64.1584/4.0],[40,10];
    ic,cfg = ic_slump(L,nel);
    status = slump(ic,cfg;);
=#