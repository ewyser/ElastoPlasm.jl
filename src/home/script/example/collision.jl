export ic_collision

"""
    ic_collision(L::Vector{Float64}, nel::Vector{Int64}; fid::String=..., kwargs...) -> String

Initializes the mesh, material points, and simulation configuration for a two-disk
elastic collision test, and exports it to a `.jld2` file (same pipeline as `slump_problem`).

# Arguments
- `L::Vector{Float64}`: Domain dimensions.
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration (see `get_solver`).

# Returns
- `String`: Path to the exported `.jld2` setup file.

# Example
```julia
jld2 = ic_collision([64.0, 16.0], [40, 10]; fid="run1")
out  = elastoplasm(jld2; workflows=[elastodynamic!])
```
"""
function ic_collision(L,nel; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d disk collision problem"
    # get solver and paths
    solver = get_solver(; dim=length(L), kwargs...)
    paths  = set_paths(fid,self.sys.out;interactive=false)
    # mesh, mpts, mat & time initial conditions
    geom   = setup_geometry(L,nel,solver)
    mesh   = setup_mesh(geom,solver)
    mat    = setup_material_constants(solver                                    )
    mpts   = setup_mpts(mesh,solver,mat ; geom = get_collision(mesh,mat,solver))
    time   = setup_time(solver     ; te = 10.0, tg = 10.0, tep = 5.0  )
    problem= MechanicalProblem(mesh,mpts,time)
    basis  = setup_basis(problem,solver)
    # plot initial velocity components
    if solver.plot.status
        @info "Plotting initial velocity components..."
        dims = solver.plot.dpi .* (mesh.prprt.L ./ mesh.prprt.L[1])
        ms   = maximum(dims) / (maximum(mesh.prprt.nel) * 0.005)
        opts = (;
            dims    = (solver.plot.dpi .* mesh.prprt.L[1], 2 * solver.plot.dpi .* mesh.prprt.L[2]),
            what    = [
                (;mpts=(name="u", cblim=(-25.0, 25.0))),
                (;mpts=(name="v", cblim=(-25.0, 25.0)))
            ],
            xlim    = (minimum(getindex.(mesh.x, 1)), maximum(getindex.(mesh.x, 1))),
            ylim    = (minimum(getindex.(mesh.x, 2)), maximum(getindex.(mesh.x, 2))),
            tit     = L" t = " * string(round(0.0, digits=1)) * " [s]",
            backend = gr(legend=true, markersize=ms, markershape=:circle, markerstrokewidth=0.75,),
            file    = joinpath(paths[:plot], "$(length(mesh.prprt.nel)-1)d_collision_velocity.png"),
        )
        get_plot_field(mpts,mesh,opts); save_plot(opts)
    end
    # export to jld2 file and return path
    return export_problem(problem,basis,solver,paths; path = paths[:dat], file = "collision_simulation")
end
