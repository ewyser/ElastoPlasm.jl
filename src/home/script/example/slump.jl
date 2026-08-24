export slump_problem 

"""
    slump_problem(L::Vector{Float64}, nel::Vector{Int64}; fid::String=..., kwargs...) -> NamedTuple, NamedTuple

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
ic, cfg = slump_problem([64.0, 16.0], [40, 10]; fid="run1")
```
"""
function slump_problem(L,nel; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d slump problem"
    # get solver and paths
    solver = get_solver(; dim=length(L), kwargs...)
    paths  = set_paths(fid,self.sys.out;interactive=false)  
    # mesh, mpts, mat & time initial conditions
    geom    = setup_geometry(L,nel,solver)
    problem = setup_problem(geom,solver,get_slump; te = 10.0, tg = 10.0, tep = 5.0)
    basis   = setup_basis(problem,solver)
    # plot initial cohesion field
    if solver.plot.status
        @info "Plotting initial cohesion & friction fields..."
        dims  = solver.plot.dpi.*(problem.mesh.prprt.L[1]./problem.mesh.prprt.L)
        ms    = dims[1]/(problem.mesh.prprt.nel[1]*2)

        what = [(;mpts=get_mpts_variable_config()[name]) for name ∈ ["coh0", "phi0"]]
        opts = (;
            dims    = solver.plot.dpi.*(problem.mesh.prprt.L./problem.mesh.prprt.L[1]),
            what    = what,
            xlim    = (minimum(getindex.(problem.mesh.x, 1)), maximum(getindex.(problem.mesh.x, 1))),
            ylim    = (minimum(getindex.(problem.mesh.x, 2)), maximum(getindex.(problem.mesh.x, 2))),
            tit     = L" t = "*string(round(0.0,digits=1))*" [s]",
            backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
            file    = joinpath(paths[:plot],"$(length(problem.mesh.prprt.nel)-1)d_std_coh0_phi0.png"),
        )
        get_plot_field(problem.mpts,problem.mesh,opts);save_plot(opts)
    end
    # export to jld2 file and return path
    return export_problem(problem,basis,solver,paths; path = paths[:dat], file = "slump_simulation")
end