export require
"""
    require(in::Symbol = :instr) -> NamedTuple

Return a reference configuration named tuple for simulation instructions or API routines.

# Arguments
- `in::Symbol`: The configuration type to load defaults for. Only `:instr` is currently supported.

# Returns
- `NamedTuple`: Default configuration values for simulation, including precision, basis, deformation framework, mapping scheme, locking mitigation, random field generator, plasticity, non-local regularization, plotting options, and performance mode.

# Example
```julia
cfg = require(:instr)
println(cfg.basis.which)  # prints the default basis type
```

# Notes
- If `in == :instr`, returns a named tuple with all default simulation options.
- Throws an error for unsupported symbols.

# Supported keys
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
function require(in::Symbol=:instr)
    if in == :instr
        return (;
            dtype = 64,
            basis = (;
                which = "bsmpm",
                how = nothing,
                ghost = false,
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
                status=true,
                freq=1.0,
                what=["T"],
                dims=(500.0,250.0),
            ),
            perf  = (;
                status=false,
            ),
        )
    else
        throw(error("UnsupportedSymbol: $(in)"))
        return nothing
    end
end