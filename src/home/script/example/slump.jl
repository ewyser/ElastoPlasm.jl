export slump, slump!, ic_slump 

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
    instr = kwargser(:instr,kwargs;dim=length(L))
    paths = set_paths(fid,info.sys.out;interactive=false)      
    # mesh & mp initial conditions
    mesh  = setup_mesh(nel,L,instr)    
    cmpr  = setup_cmpr(mesh.dim,instr)                       
    mp    = setup_mps(mesh,cmpr;define=inislump(mesh,cmpr,instr))
    # time parameters
    te,tg = 10.0, 10.0
    tep   = 5.0
    time  = (; t = [0.0,te+tep], te = te, tg = if tg > te te else tg end, tep = tep,)
    # plot initial cohesion field
    plotcoh(mp,cmpr,paths)   
    # display summary
    @info ic_log(mesh,mp,time)
    return (;mesh,mp,cmpr,time),(;instr,paths)
end

"""
    slump(ic::NamedTuple, cfg::NamedTuple; mutate::Bool=false) -> Bool

Runs the explicit solution workflow for the slump problem, including simulation and postprocessing.

# Arguments
- `ic::NamedTuple`: Initial mesh, material point, compression, and time configuration.
- `cfg::NamedTuple`: Simulation instructions and output paths.

# Returns
- `Bool`: Returns `true` if the simulation and postprocessing complete successfully.

# Example
```julia
success = slump(ic, cfg)
```
"""
function slump(ic::NamedTuple,cfg::NamedTuple; load::String="elastodynamic")
    @info "Explicit solution to slump problem";config_plot()                                           
    # forward-euler explicit workflow
    out = elastoplasm(deepcopy(ic),deepcopy(cfg); mode = load)
    # return success
    return success=true
end
function slump!(ic::NamedTuple,cfg::NamedTuple; load::String="elastodynamic")
    @info "Explicit solution to slump problem";config_plot()                                           
    # forward-euler explicit workflow
    out = elastoplasm(ic          ,cfg         ; mode = load)
    # return success
    return out = (; out...,success=true,)
end
#=
    L,nel  = [64.1584,64.1584/4.0],[40,10];
    ic,cfg = ic_slump(L,nel);
    status = slump(ic,cfg;);
=#