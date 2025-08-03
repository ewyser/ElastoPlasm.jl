export slump,slump!,ic_slump 

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
    cmpr  = setup_cmpr(mesh,instr)                       
    mp    = setup_mps(mesh,cmpr;define=geom_slump(mesh,cmpr,instr))
    # time parameters
    te,tg = 10.0, 10.0
    tep   = 5.0
    time  = (; t = [0.0,te+tep], te = te, tg = if tg > te te else tg end, tep = tep,)
    # plot initial cohesion field
    ms = 0.4*instr[:plot][:dims][1]/mesh.nel[1]
    opts = (;
        dims    = instr[:plot][:dims],
        what    = ["coh0","phi0"],
        backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
        tit     = L" t = "*string(round(0.0,digits=1))*" [s]",
        file    = joinpath(paths[:plot],"$(mesh.dim)d_coh0_phi0.png"),
    )
    get_plot_field(mp,mesh,opts);save_plot(opts)
    # display summary
    @info ic_log(mesh,mp,time)
    return (;mesh,mp,cmpr,time),(;instr,paths)
end

"""
    slump(ic::NamedTuple, cfg::NamedTuple; workflow::String="elastodynamic") -> NamedTuple

Runs the explicit solution workflow for the slump problem using a deep copy of the initial conditions and configuration.
This function is suitable for workflows where you do not want to mutate the input data.

# Arguments
- `ic::NamedTuple`: Initial mesh, material point, compression, and time configuration.
- `cfg::NamedTuple`: Simulation instructions and output paths.
- `workflow::String`: (Optional) Workflow mode, default is "elastodynamic".

# Returns
- `NamedTuple`: Simulation output with all fields from `elastoplasm` and an added `success=true` field.

# Example
```julia
result = slump(ic, cfg)
if result.success
    println("Simulation completed successfully!")
end
```
"""
function slump(ic::NamedTuple,cfg::NamedTuple; workflow::String="elastodynamic")
    @info "Explicit solution to slump problem"; config_plot()
    # forward-euler explicit workflow
    out = elastoplasm(deepcopy(ic), deepcopy(cfg); mode = workflow)
    # return output with success flag
    return out = (; out..., success=true,)
end

"""
    slump!(ic::NamedTuple, cfg::NamedTuple; workflow::String="elastodynamic") -> NamedTuple

Runs the explicit solution workflow for the slump problem, mutating the input initial conditions and configuration.
Use this when you want changes to `ic` and `cfg` to persist after the simulation.

# Arguments
- `ic::NamedTuple`: Initial mesh, material point, compression, and time configuration.
- `cfg::NamedTuple`: Simulation instructions and output paths.
- `workflow::String`: (Optional) Workflow mode, default is "elastodynamic".

# Returns
- `NamedTuple`: Simulation output with all fields from `elastoplasm` and an added `success=true` field.

# Example
```julia
result = slump!(ic, cfg)
if result.success
    println("Simulation completed successfully!")
end
```
"""
function slump!(ic::NamedTuple,cfg::NamedTuple; workflow::String="elastodynamic")
    @info "Explicit solution to slump problem"; config_plot()
    # forward-euler explicit workflow
    out = elastoplasm(ic, cfg; mode = workflow)
    # return output with success flag
    return out = (; out..., success=true,)
end