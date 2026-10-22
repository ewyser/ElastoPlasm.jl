export plot_initial_fields

"""
    plot_initial_fields(problem, folder)

Plot the initial cohesion (and, for Drucker–Prager, friction angle) of the material points into
`folder`, as `setup_simulation` does when `solver.plot.status` is on.
"""
function plot_initial_fields(problem::Problem, folder)
    solver = problem.solver
    @info "Plotting initial cohesion & friction fields..."
    dims  = solver.plot.dpi.*(problem.mesh.prprt.L[1]./problem.mesh.prprt.L)
    ms    = dims[1]/(problem.mesh.prprt.nel[1]*2)

    # phi0 (friction angle) only exists on DruckerPrager's constitutive constants —
    # VonMises (material.plastic=="VM") has no friction angle at all (J2 yield has
    # no pressure dependence, see VonMises's docstring), so plotting it there would
    # either crash (no .ϕ₀ field) or fake a value that doesn't mean anything.
    what  = solver.material.plastic == "DP" ? ["c0", "phi0"] : ["c0"]
    opts = (;
        dims    = solver.plot.dpi.*(problem.mesh.prprt.L./problem.mesh.prprt.L[1]),
        what    = what,
        xlim    = (minimum(getindex.(problem.mesh.x, 1)), maximum(getindex.(problem.mesh.x, 1))),
        ylim    = (minimum(getindex.(problem.mesh.x, 2)), maximum(getindex.(problem.mesh.x, 2))),
        tit     = L" t = "*string(round(0.0,digits=1))*" [s]",
        backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
        file    = joinpath(folder,"$(length(problem.mesh.prprt.nel)-1)d_std_coh0_phi0.png"),
    )
    get_plot_field(problem.mpts,problem.mesh,opts);save_plot(opts)
end
