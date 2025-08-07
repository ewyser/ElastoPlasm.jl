export require
"""
    require(in::Symbol=:instr) -> Dict

Returns a reference configuration dictionary for simulation instructions or API routines.

# Arguments
- `in::Symbol`: The configuration type to load defaults for. Only `:instr` is supported.

# Returns
- `Dict`: A dictionary containing default configuration values for simulation, including precision, basis, deformation framework, mapping scheme, locking mitigation, random field generator, plasticity, non-local regularization, plotting options, and performance mode.

# Behavior
- If `in == :instr`, returns a dictionary with all default simulation options.
- Throws an error for unsupported symbols.

# Example
```julia
cfg = require(:instr)
```

# Supported keys and their purpose
- `:dtype`   — Arithmetic precision
- `:basis`   — Shape function type
- `:fwrk`    — Deformation framework
- `:trsfr`   — Mapping scheme
- `:vollock` — Volumetric locking mitigation
- `:grf`     — Gaussian Random Field generator
- `:plast`   — Plasticity onset and flow law
- `:nonloc`  — Non-local regularization
- `:plot`    — Plotting options
- `:perf`    — Performance mode
"""
function require(in::Symbol=:instr)
    if in == :instr
        instr = Dict(
            :dtype => 64,
            :basis => (;
                        which="bsmpm",
                        how=nothing,
                        ghost=false,
            ),
            :fwrk  => (;
                        deform = "finite",
                        trsfr = "musl",
                        locking = true,
            ),
            :bcs   => (;
                        dirichlet = :roller,
            ),
            :grf   => (;
                        status = false,
                        covariance = "gaussian",
                        param = (; Iₓ= [2.5,2.5,2.5], Nₕ = 5000, kₘ = 100,),
            ),
            :plast => (;
                        status = false,
                        constitutive = "DP",
            ),
            :nonloc=> (;
                        status=true,
                        ls=0.5,
            ),
            :plot  => (;
                        status=true,
                        freq=1.0,
                        what=["epII"],
                        dims=(500.0,250.0),
            ),
            :perf    => false,
        )
        return instr
    else
        throw(error("UnsupportedSymbol: $(in)"))
        return nothing
    end
end