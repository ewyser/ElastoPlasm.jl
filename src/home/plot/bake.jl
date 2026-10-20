"""
    bake(mpts, mesh, t, solver; title="t = … [s]")

Plot and display the fields of `solver.plot.what` when `solver.plot.status`; each panel is titled
"<field>, at <title>", by default the time `t` (a workflow without physical time, such as the
load steps of `elastoquasistatic!`, passes its own `title`).
"""
@views function bake(mpts::Point{T1,T2},mesh::Mesh{T1,T2},t::T2,solver::S; title=L" t = "*string(round(t,digits=1))*" [s]") where {T1,T2,S<:AbstractSolver{T1,T2}}
    if solver.plot.status
        opts = (;
            dims    = solver.plot.dpi.*(mesh.prprt.L./mesh.prprt.L[1]),
            what    = solver.plot.what,
            xlim    = (minimum(getindex.(mesh.x, 1)), maximum(getindex.(mesh.x, 1))),
            ylim    = (minimum(getindex.(mesh.x, 2)), maximum(getindex.(mesh.x, 2))),
            tit     = title,
            backend = gr(legend=true,markershape=:circle,markerstrokewidth=0.75,),
        )
        return get_plot_field(mpts,mesh,opts)
    else
        nothing
    end
    return nothing
end