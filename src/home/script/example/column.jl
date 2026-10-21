export column_case

"""
    column_case(nel, ν, E, ρ0, l0; folder="column", kwargs...) -> Case

An elastic self-weight column of height `l0` in a domain `1.25 l0` high, one element wide (2D, `nel =
[1, nz]`, or 3D, `[1, 1, nz]`), with Poisson's ratio `ν`, Young's modulus `E` and density `ρ0`: the
convergence test of `test/testset/test_column.jl` against the geostatic stress. Gravity is ramped over
`tg = ⌈40 · 2 l0 / c⌉` (`c` the elastic wave speed), which is also the end time. `kwargs` are solver
options. Not the granular `collapse_case`.

```julia
jld2 = setup_simulation(column_case([1, 20], 0.0, 1.0e4, 80.0, 50.0))
execute!(jld2; workflows=[elastodynamic!])
```
"""
function column_case(nel, ν, E, ρ0, l0; folder::String="column", kwargs...)
    D  = length(nel)
    ly = 1.25*l0
    lx = ly/nel[end]
    L  = D == 2 ? [lx, ly] : [lx, lx, ly]
    K, G = get_elastic_moduli(Float64(E), Float64(ν))
    tg = ceil((1.0/sqrt((K + 4.0/3.0*G)/Float64(ρ0)))*(2.0*l0)*40.0)
    return Case(; L, nel, ni=2, zones=[Zone("column", Block(zeros(D), [L[1:end-1]; l0]))],
                  material=(; E, nu=ν, rho0=ρ0), solver=(; kwargs...), run=(; folder, te=tg, tg=tg, tep=0.0))
end
