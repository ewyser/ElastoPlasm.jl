export get_default

"""
    get_default() -> NamedTuple

Return default configuration values for the solver.

# Returns
- `NamedTuple`: Default configuration values for simulation, including precision, basis, deformation framework, mapping scheme, locking mitigation, random field generator, plasticity, non-local regularization, plotting options, and performance mode.

# Example
```julia
cfg = get_default()
println(cfg.basis.which)  # prints the default basis type
```

# Configuration Keys
- `:dtype`    — Arithmetic precision (e.g., 64 for Float64)
- `:basis`    — Shape function type/options, and the P2G/G2P transfer scheme and its blend knob
- `:strain`   — Strain formulation (finite/infinitesimal)
- `:stab`     — Numerical stabilization (F-bar locking correction, damping, MUSL reprojection)
- `:bcs`      — Boundary condition settings
- `:grf`     — Gaussian Random Field generator options
- `:plast`   — Plasticity onset and flow law
- `:nonloc`  — Non-local regularization options
- `:plot`    — Plotting options
- `:perf`    — Performance mode options
"""
function get_default()
    default = (;
        solution = "explicit",
        dtype = (;
            T0 = (Int64,Float64),
            bits = Int64(64),
            precision = "64-bit precision (or double precision)"
        ),
        basis = (;
            which = "bsmpm",
            how = nothing,
            trsfr = "std",
            C_pf = 1.0,
        ),
        strain   = (;
            deform = "finite",
        ),
        stab     = (;
            locking = true,
            damping = 0.1,
            musl = true,
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
            what   = [(;mpts=get_mpts_variable_config()["epII"]),],
        ),
        perf  = (;
            status=false,
        ),
        backend = (;
            select="host",
            distributed=false
        ),
    )
    return default
end
