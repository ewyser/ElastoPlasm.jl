export ic_slump 

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
function ic_slump(L,nel; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d slump problem"
    # init & kwargs
    instr = kwargser(kwargs; dim=length(L))
    paths = set_paths(fid,info.sys.out;interactive=false)  
    # mesh, mpts, cmpr & time initial conditions
    geom  = setup_geometry(L,nel,instr)
    mesh  = setup_mesh(geom,instr)
    cmpr  = setup_cmpr(mesh                                         )                       
    mpts  = setup_mpts(mesh,instr,cmpr ; geom = get_slump(mesh,cmpr,instr))
    time  = setup_time(instr     ; te = 10.0, tg = 10.0, tep = 5.0  ) 
    # plot initial cohesion field
    if instr.plot.status
        @info "Plotting initial cohesion & friction fields..."
        dims  = instr.plot.dpi.*(mesh.prprt.L[1]./mesh.prprt.L)
        ms    = dims[1]/(mesh.prprt.nel[1]*2)

        what = [(;mpts=get_mpts_variable_config()[name]) for name ∈ ["coh0", "phi0"]]
        opts = (;
            dims    = instr.plot.dpi.*(mesh.prprt.L./mesh.prprt.L[1]),
            what    = what,
            xlim    = (minimum(mesh.x[1,:]),maximum(mesh.x[1,:])),
            ylim    = (minimum(mesh.x[2,:]),maximum(mesh.x[2,:])),
            tit     = L" t = "*string(round(0.0,digits=1))*" [s]",
            backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
            file    = joinpath(paths[:plot],"$(mesh.prprt.dim)d_std_coh0_phi0.png"),
        )
        get_plot_field(mpts,mesh,opts);save_plot(opts)
    end
    # display summary
    @info ic_log(mesh,mpts,time,instr)
    misc = (;
        prefix = "$(mesh.prprt.dim)d_$(instr.fwrk.trsfr)"
    )
    # export to jld2 file and return path
    return export_setup(mesh,mpts,cmpr,time,instr,paths,misc; path = paths[:dat], file = "slump_simulation")
end