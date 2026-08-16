export ic_collision 

"""
    ic_collision(L::Vector{Float64}, nel::Vector{Int64}; fid::String=..., kwargs...) -> NamedTuple, NamedTuple

Initializes the mesh, material points, and simulation configuration for a disk collision test.

# Arguments
- `L::Vector{Float64}`: Domain dimensions.
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration.

# Returns
- `(ic, cfg)`: Two named tuples containing mesh/material point/compression info (`ic`) and instructions/paths (`cfg`).

# Example
```julia
ic, cfg = ic_collision([64.0, 16.0], [40, 10]; fid="run1")
```
"""
function ic_collision(L, nel; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d disk collision problem"
    # init & kwargs
    instr = kwargser(kwargs; dim=length(L))
    paths = set_paths(fid, self.sys.out; interactive=false)  
    # mesh, mpts, cmpr & time initial conditions
    mesh  = setup_mesh(instr     ; geom = get_geom(nel, L, instr)     )
    cmpr  = setup_cmpr(mesh                                           )                       
    mpts  = setup_mpts(mesh, instr, cmpr ; geom = get_collision(mesh, cmpr, instr))
    time  = setup_time(instr     ; te = 10.0, tg = 10.0, tep = 5.0   ) 
    # plot initial velocity components
    if instr.plot.status
        @info "Plotting initial velocity components..."
        dims  = instr.plot.dpi .* (mesh.prprt.L ./ mesh.prprt.L[1])
        ms    = maximum(dims) / (maximum(mesh.prprt.nel) * 0.005)
        # Plot both velocity components in a single figure
        opts = (;
            dims    = (instr.plot.dpi .* mesh.prprt.L[1], 2 * instr.plot.dpi .* mesh.prprt.L[2]),
            what    = [
                (;mpts=(name="u", cblim=(-25.0, 25.0))),
                (;mpts=(name="v", cblim=(-25.0, 25.0)))
            ],
            xlim    = (mesh.prprt.xB[1, 1], mesh.prprt.xB[1, 2]),
            ylim    = (mesh.prprt.xB[2, 1], mesh.prprt.xB[2, 2]),
            tit     = L" t = " * string(round(0.0, digits=1)) * " [s]",
            backend = gr(legend=true, markersize=ms, markershape=:circle, markerstrokewidth=0.75,),
            file    = joinpath(paths[:plot], "$((length(mesh.prprt.nel)-1))d_collision_velocity.png"),
        )
        get_plot_field(mpts, mesh, opts); save_plot(opts)
    end
    # display summary
    @info ic_log(mesh, mpts, time, instr)
    misc = (;
        prefix = "$((length(mesh.prprt.nel)-1))d_$(instr.fwrk.trsfr)_collision"
    )
    # export to jld2 file and return path
    return export_setup(mesh, mpts, cmpr, time, instr, paths, misc; path = paths[:dat], file = "collision_simulation")
end
