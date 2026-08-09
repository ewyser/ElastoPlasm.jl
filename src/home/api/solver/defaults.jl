export get_default

"""
    get_default(::Type{Instruction}) -> NamedTuple

Return default configuration values for the specified type.

# Arguments
- `::Type{Instruction}`: The type to get defaults for

# Returns
- `NamedTuple`: Default configuration values for simulation, including precision, basis, deformation framework, mapping scheme, locking mitigation, random field generator, plasticity, non-local regularization, plotting options, and performance mode.

# Example
```julia
cfg = get_default(Instruction)
println(cfg.basis.which)  # prints the default basis type
```

# Configuration Keys
- `:dtype`   — Arithmetic precision (e.g., 64 for Float64)
- `:basis`   — Shape function type and options
- `:fwrk`    — Deformation framework and transfer scheme
- `:bcs`     — Boundary condition settings
- `:grf`     — Gaussian Random Field generator options
- `:plast`   — Plasticity onset and flow law
- `:nonloc`  — Non-local regularization options
- `:plot`    — Plotting options
- `:perf`    — Performance mode options
"""
function get_default()
    solver  = ExplicitSolver
    default = (;
        dtype = (;
            T0=(Int64,Float64),
            bits=Int64(64),
            precision="64-bit precision (or double precision)"
        ),
        basis = (;
            which = "bsmpm",
            how = nothing,
            ghost = 0,
        ),
        fwrk  = (;
            deform = "finite",
            trsfr = "std",
            C_pf = 1.0, 
            musl = true,
            locking = true,
            damping = 0.1
        ),
        bcs   = (;
            dirichlet = [
                :roller :roller;
                :roller :roller;
                :roller :roller], # for 2d, this translates to [lower_x upper_x;lower_y upper_y]
        ),
        grf   = (;
            status = false,
            covariance = "gaussian",
            param = (; 
                Iₓ= [2.5,2.5,2.5], 
                Nₕ = 5000, 
                kₘ = 100,
            ),
        ),
        plast = (;
            status = false,
            constitutive = "DP",
        ),
        nonloc = (;
            status=true,
            ls=0.5,
        ),
        plot  = (;
            status = true,
            freq   = 1.0,
            dpi    = 500,
            what   = [(;mpts=get_mpts_variable_config()["P"]),],
        ),
        perf  = (;
            status=false,
        ),
        backend = (;
            select="host",
            distributed=false
        ),
    )
    return solver, default
end
