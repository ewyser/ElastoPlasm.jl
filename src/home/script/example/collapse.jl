export collapse_case

"""
    collapse_case(nel; L=[0.92, 0.22], w0=0.2, h0=0.1, ρ0=2650.0, E, ν=0.3, ϕ=deg2rad(19.8), c0=0.0, folder="collapse", workflow, kwargs...) -> Case

A dry granular block collapse: a block `[0, w0] × [0, h0]` at the domain's lower-left corner,
released under gravity, Drucker–Prager, on a `:fixed` (no-slip) base, as MaterialPointSolver.jl's
`2d_collapse.jl` reference (`LandslideSIM/Archive_MaterialPointSolver.jl_paper`): domain, block and
material (`ρ0`, `E` from `Ks = 7e5`, `ν`, `ϕ`, cohesionless `c0`, `cr = c0`) match it; `nel` is the
argument meant to vary. `workflow` is the run's (`RUN.workflow` by default); `kwargs` are solver options over
those defaults. Not to be confused with
`column_case`, the elastic self-weight column.

```julia
jld2 = setup_simulation(collapse_case([92, 22]))
execute!(jld2; workflows=[elastodynamic!, elastoplastic!])
```
"""
function collapse_case(nel; L=[0.92, 0.22], w0=0.2, h0=0.1, ρ0=2650.0, E=7.0e5*3.0*(1.0-2.0*0.3), ν=0.3, ϕ=deg2rad(19.8), c0=0.0,
                       folder::String="collapse", workflow::String=RUN.workflow, kwargs...)
    solver = merge((; material=(; plastic="DP", elastic="hencky"), bcs=(; dirichlet=[:roller :roller; :fixed :roller])), (; kwargs...))
    return Case(; L, nel, ni=2, zones=[Zone("block", Block([0.0, 0.0], [w0, h0]))],
                  material=(; E, nu=ν, rho0=ρ0, c0, cr=c0, phi=ϕ), solver, run=(; folder, workflow, te=10.0, tg=10.0, tep=5.0))
end
