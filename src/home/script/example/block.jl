 export block_problem

"""
    block_problem(nel=[16, 9]; kwargs...) -> String

The 8 m × 8 m self-weight block of Sutcliffe et al. (2026, IJNME 127:e70384) §4.3, Fig 10a, built
with `collapse_problem`: `ρ = 1000 kg/m³` (`ρ0 = 1000/0.9`, since the mass is `(1-n₀)ρ₀Ω` with
`n₀ = 0.1`), `E = 1 MPa`, `ν = 0.3`, perfectly plastic von Mises with `c0 = ρ_y = 20 kPa` (the
Lode-radius reading, see `architecture.md`; `c0 = √2·20e3` for `√J₂ = ρ_y` as Eq 49 prints it), on a
16 m × 9 m grid (`h = 1 m` with the default `nel`), 2×2 particles per cell, smooth base, GIMP.
Every value is a keyword override (`L`, `w0`, `h0`, `ρ0`, `E`, `ν`, `c0`, or any `get_solver`
section). Solve it with dynamic relaxation and the paper's gravity, elastic or plastic:

```julia
jld2 = block_problem()
elastoplasm!(jld2; workflows=[(a...) -> elastoquasistatic!(a...; lstps=10, G=10.0, plastic=false)])
```
"""
function block_problem(nel=[2*16, 2*9]; kwargs...)
    defaults = (;
        L = [16.0, 9.0], w0 = 8.0, h0 = 8.0,
        ρ0 = 1000.0/0.9, E = 1.0e6, ν = 0.3, c0 = 20.0e3,
        fid = "block",
        material = (; plastic="VM", elastic="hencky"),
        bcs      = (; dirichlet=[:roller :roller; :roller :roller]),
        basis    = (; get_default().basis..., which="gimpm"),
    )
    return collapse_problem(nel; merge(defaults, kwargs)...)
end
