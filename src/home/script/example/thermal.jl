export thermal_case

"""
    thermal_case(L, nel; folder="thermal", kwargs...) -> Case

Heat conduction: the whole domain filled with material at the initial temperature, every side
`:fixed`, the thermal phase on, plotting the temperature. `kwargs` are solver options over those
defaults.

```julia
jld2 = setup_simulation(thermal_case([16.0, 16.0], [80, 80]))
execute!(jld2; workflows=[thermodynamic!])
```
"""
function thermal_case(L, nel; folder::String="thermal", kwargs...)
    D = length(L)
    solver = merge((;
            basis = (; get_default().basis..., which = "bsmpm"),
            bcs   = (; dirichlet = fill(:fixed, D, 2)),
            plot  = (; status = true, freq = 1.0, dpi = 300, what = [(;mpts=get_mpts_variable_config()["T"])]),
        ), (; kwargs...))
    return Case(; L, nel, ni=2, zones=[Zone("material", Block(zeros(D), L))], thermal=true, solver,
                  run=(; folder, te=5.0, tg=5.0, tep=0.0))
end
