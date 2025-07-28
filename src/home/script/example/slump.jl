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
    @info """
    Summary: 
    - elements: $(mesh.nel[end])
    - material points: $(mp.nmp) 
    - simulation time ∈ $(time.t) s:
        - gravity ramp-up: $(time.tg ) s    
        - elastodynamic  : $(time.te ) s
        - elastoplastic  : $(time.tep) s
    """
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
    # unpack mesh, mp, cmpr, instr, paths
    mesh,mp,cmpr = deepcopy(ic[:mesh]  ), deepcopy(ic[:mp]    ), deepcopy(ic[:cmpr])
    instr,paths  = deepcopy(cfg[:instr]), deepcopy(cfg[:paths])
    time         = deepcopy(ic[:time]  )                                             
    # action
    elastoplasm(mp,mesh,cmpr,time,paths,instr; load = load)
    # return success
    return success=true
end
function slump!(ic::NamedTuple,cfg::NamedTuple; load::String="elastodynamic")
    @info "Explicit solution to slump problem";config_plot()
    # unpack mesh, mp, cmpr, instr, paths
    mesh,mp,cmpr = ic[:mesh]  , ic[:mp]    , (ic[:cmpr])
    instr,paths  = cfg[:instr], cfg[:paths]
    time         = ic[:time]                                                
    # action
    elastoplasm(mp,mesh,cmpr,time,paths,instr; load = load)
    # return success
    return success=true
end
#=
    L,nel  = [64.1584,64.1584/4.0],[40,10];
    ic,cfg = ic_slump(L,nel);
    status = slump(ic,cfg;);
=#