"""
    savlot(mpts::Point{T1,T2,E,R}, mesh::Mesh{T1,T2,D}, t::T2, instr::NamedTuple) where {D}

Plot and display simulation fields at the current time step, if plotting is enabled in the instruction dictionary.

# Arguments
- `mpts::Point{T1,T2,E,R}`: Material point data structure.
- `mesh::Mesh{T1,T2,D}`: Mesh data structure.
- `t::T2`: Current simulation time.
- `instr::NamedTuple`: Simulation instruction dictionary, must contain plotting options in `:plot`.

# Returns
- `Nothing`. Displays the plot(s) if plotting is enabled, otherwise returns `nothing`.

# Example
```julia
savlot(mpts, mesh, t, instr)
```
"""
@views function savlot(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},t::T2,instr::Instruction{T1,T2,D}) where {T1,T2,E,R,D}
    if instr.plot.status
        dims = instr.plot.dpi.*(mesh.prprt.L[1]./mesh.prprt.L)
        ms   = dims[1]/(mesh.prprt.nel[1]*2)
        opts = (;
            dims    = instr.plot.dpi.*(mesh.prprt.L./mesh.prprt.L[1]), 
            what    = instr.plot.what,
            xlim    = (minimum(mesh.x[1,:]),maximum(mesh.x[1,:])), 
            ylim    = (minimum(mesh.x[2,:]),maximum(mesh.x[2,:])),
            cblim   = instr.plot.cblim,
            tit     = L" t = "*string(round(t,digits=1))*" [s]",
            backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
        )
        return get_plot_field(mpts,mesh,opts)
    else
        nothing
    end
    return nothing
end