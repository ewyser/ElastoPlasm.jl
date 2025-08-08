"""
    savlot(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, t::T2, instr::NamedTuple)

Plot and display simulation fields at the current time step, if plotting is enabled in the instruction dictionary.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `t::T2`: Current simulation time.
- `instr::NamedTuple`: Simulation instruction dictionary, must contain plotting options in `:plot`.

# Returns
- Displays the plot(s) if plotting is enabled, otherwise returns `nothing`.

# Example
```julia
savlot(mpts, mesh, t, instr)
```
"""
@views function savlot(mpts::Point{T1,T2},mesh::Mesh{T1,T2},t::T2,instr::NamedTuple) where {T1,T2}
    if instr[:plot][:status]
        ms = 0.4*instr[:plot][:dims][1]/mesh.nel[1]
        opts = (;
            dims  = instr[:plot][:dims],
            what = instr[:plot][:what],
            backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
            tit     = L" t = "*string(round(t,digits=1))*" [s]",
        )
        return get_plot_field(mpts,mesh,opts)
    else
        nothing
    end
    return nothing
end