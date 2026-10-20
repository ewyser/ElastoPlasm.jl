export slump_problem, slump_case

"""
    slump_case(L, nel; lz=12.8, ni=2, fid="slump", kwargs...) -> Case
    slump_problem(L, nel; kwargs...) -> String

The slump: a slope `z = 0.625 Lx − 1.25 x` cut into a layer of height `lz` on a base layer
`z < 0.15 lz` that spans the domain; below `z = 0.3 lz` the material is weaker (friction angle 7.5°)
and three times warmer. `L`, `nel` are 2D (x, z) or 3D (x, y, z, extruded along y); the material
fills the sub-cells wholly below `lz` (the slope polygon's top is on the sub-cell lattice). `kwargs`
are solver options (`get_solver`), e.g. `cli()...`; the random field `grf` perturbs the cohesion.
`slump_problem` builds it (`setup_simulation`) and returns the `.jld2` path.

```julia
jld2 = slump_problem([64.1584, 16.0396], [40, 10]; cli()...)
elastoplasm!(jld2; workflows=[elastodynamic!, elastoplastic!])
```
"""
function slump_case(L, nel; lz=12.80, ni::Integer=2, fid::String="slump", kwargs...)
    Lx, sub = L[1], L[end]/(nel[end]*ni)                      # sub-cell height
    wl   = 0.15*lz
    top  = min(sub*floor(lz/sub), 0.625*Lx)                    # top of the whole sub-cells below lz (or of the slope)
    xs(z) = (0.625*Lx - z)/1.25                                # the slope line
    slope = Polygon([[0.0, 0.0], [Lx, 0.0], [Lx, wl], [xs(wl), wl], [xs(top), top], [0.0, top]])
    base  = Polygon([[0.0, 0.0], [Lx, 0.0], [Lx, wl], [xs(wl), wl], [xs(2wl), 2wl], [0.0, 2wl]])
    zones = [Zone("slope", slope), Zone("base", base; phi=7.5*π/180, T=3.0*MATERIAL.T)]
    return Case(; L, nel, ni, zones, solver=(; kwargs...), run=(; fid, te=10.0, tg=10.0, tep=5.0))
end
slump_problem(L, nel; kwargs...) = setup_simulation(slump_case(L, nel; kwargs...))
