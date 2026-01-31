export get_option

"""
    get_option(::Type{Instruction}) -> NamedTuple

Return all possible values for each configuration field as vectors.

# Arguments
- `::Type{Instruction}`: The type to get options for

# Returns
- `NamedTuple`: Structure matching `get_default(Instruction)` but with vectors of possible values instead of single values.

# Example
```julia
opts = get_option(Instruction)
println(opts.basis.which)  # ["bsmpm", "gimpm", "smpm"]
println(opts.fwrk.deform)  # ["finite", "infinitesimal"]
```

# Notes
- Useful for validation, UI generation, and documentation
- Fields with continuous values (like damping) show reasonable ranges
- Boolean fields are represented as `[true, false]`
"""
function get_option(::Type{Instruction})
    return (
            dtype = ("Select arithmetic precision",[32, 64]),
            basis = (
                which = ("Select basis type",["bsmpm", "gimpm", "smpm"]),
                how = ("Select material point domain update",[nothing]),
                ghost = ("Add ghost nodes ?",[true, false]),
            ),
            fwrk  = (
                deform = ("Select the deformation framework",["finite", "infinitesimal"]),
                trsfr = ("Select the mapping scheme",["std", "tpic", "apic"]),
                C_pf = ("Select picflip ratio",[1.0, 0.99, 0.95]),
                musl = ("Enable musl update",[true, false]),
                locking = ("Enable volumetric locking mitigation",[true, false]),
                damping = ("Select damping coefficient",[0.0, 0.1, 0.2, 0.4])
            ),
            grf   = (
                status = ("Enable Gaussian Random Field ",[true, false]),
                covariance = ("Select covariance function",["gaussian", "exponential"]),
                param = ( 
                    Iₓ = ("Select correlation length",[[1.0,1.0,1.0], [2.5,2.5,2.5], [5.0,5.0,5.0]]), 
                    Nₕ = ("Select number of points",[1000, 5000, 10000]), 
                    kₘ = ("Select maximum wavenumber",[50, 100, 200]),
                ),
            ),
            plast = (
                status = ("Enable plasticity",[true, false]),
                constitutive = ("Select constitutive model",["DP", "MC", "VM"]),
            ),
            nonloc = (
                status = ("Enable nonlocal effects",[true, false]),
                ls = ("Select length scale",[0.125 ,0.25, 0.5]),
            ),
            plot  = (
                status = ("Enable plotting",[true, false]),
                freq   = ("Select plot frequency",[0.1, 0.5, 1.0, 5.0]),
                dpi    = ("Select plot resolution",[100, 300, 500, 1000]),
                what   = ("Select plot variable(s)", get_variable_plot_options()),
            ),
            perf  = (
                status = ("Enable performance monitoring",[true, false]),
            ),
            backend = (
                select = ("Select backend",["host", "cuda", "rocm"]),
                distributed = ("Enable distributed computing",[true, false])
            ),
        )
end
