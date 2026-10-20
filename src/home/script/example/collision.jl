export ic_collision, collision_case

"""
    collision_case(L, nel; r=0.1, v=10.5, fid="collision", kwargs...) -> Case
    ic_collision(L, nel; kwargs...) -> String

Two elastic disks (2D) or spheres (3D) of radius `r·Lz`, at a quarter and three quarters of the
domain (3D: 0.4 and 0.6 along x and z, mid-width along y), moving towards each other at speed `v`
along the line between their centres: an SPH benchmark for momentum conservation. `kwargs` are
solver options.

```julia
jld2 = ic_collision([1.0, 1.0], [20, 20])
elastoplasm(jld2; workflows=[elastodynamic!])
```
"""
function collision_case(L, nel; r=0.1, v=10.5, fid::String="collision", kwargs...)
    D, lz = length(L), L[end]
    c0, c1 = D == 2 ? ([0.25*L[1], 0.25*lz], [0.75*L[1], 0.75*lz]) : ([0.4*L[1], 0.5*L[2], 0.4*lz], [0.6*L[1], 0.5*L[2], 0.6*lz])
    d    = c1 .- c0
    n    = sqrt(sum(d[i]*d[i] for i ∈ eachindex(d)))
    zones = [Zone("disk 1", Ball(c0, r*lz); v=v .* (d ./ n)), Zone("disk 2", Ball(c1, r*lz); v=-v .* (d ./ n))]
    return Case(; L, nel, ni=2, zones, solver=(; kwargs...), run=(; fid, te=10.0, tg=10.0, tep=5.0))
end
ic_collision(L, nel; kwargs...) = setup_simulation(collision_case(L, nel; kwargs...))
