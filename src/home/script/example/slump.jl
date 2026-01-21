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
function ic_slump(L,nel; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d slump problem"
    # init & kwargs
    instr = kwargser(kwargs, Instruction; dim=length(L))
    println(instr)
    paths = set_paths(fid,info.sys.out;interactive=false)  
    # mesh, mpts, cmpr & time initial conditions
    mesh  = setup_mesh(instr     ; geom = get_geom(nel,L,instr)     )
    cmpr  = setup_cmpr(mesh                                         )                       

println(typeof(mesh))
println(typeof(instr))
println(typeof(cmpr))

    mpts  = setup_mpts(mesh,instr,cmpr ; geom = get_slump(mesh,cmpr,instr))
    time  = setup_time(instr     ; te = 10.0, tg = 10.0, tep = 5.0  ) 
    # plot initial cohesion field
    dims  = instr.plot.dpi.*(mesh.prprt.L[1]./mesh.prprt.L)
    ms    = mesh.prprt.nel[1]*dims[1]
    opts = (;
        dims    = instr.plot.dpi.*(mesh.prprt.L./mesh.prprt.L[1]), 
        what    = [("mpts","coh0"),("mpts","phi0")],
        cblim   = [1e-3.*(cmpr.c0-cmpr.c0/2,cmpr.c0+cmpr.c0/2),(cmpr.ϕr,cmpr.ϕ0),],
        xlim    = (minimum(mesh.x[1,:]),maximum(mesh.x[1,:])), 
        ylim    = (minimum(mesh.x[2,:]),maximum(mesh.x[2,:])),        
        backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
        tit     = L" t = "*string(round(0.0,digits=1))*" [s]",
        file    = joinpath(paths[:plot],"$(mesh.prprt.dim)d_coh0_phi0.png"),
    )
    get_plot_field(mpts,mesh,opts);save_plot(opts)
    # display summary
    @info ic_log(mesh,mpts,time,instr)
    misc = (;
        prefix = "$(mesh.prprt.dim)d_$(instr.fwrk.trsfr)"
    )
    # export to jld2 file and return path
    return export_setup(mesh,mpts,cmpr,time,instr,paths,misc; path = paths[:dat], file = "slump_simulation")
end