export collapse,collapse!,ic_collapse

"""
    ic_collapse(nel::Vector{Int64}, ν, E, ρ0, l0; fid::String=..., kwargs...) -> NamedTuple, NamedTuple

Initializes the mesh, material points, constitutive model, and simulation configuration for a collapse test.

# Arguments
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `ν`: Poisson's ratio.
- `E`: Young's modulus.
- `ρ0`: Initial density.
- `l0`: Characteristic length.
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration.

# Returns
- `(ic, cfg)`: Two named tuples containing mesh/material/constitutive/time data structure (`ic`) and instructions/paths (`cfg`).

# Example
```julia
ic, cfg = ic_collapse([5, 10], 0.0, 1.0e4, 80.0, 10.0; plot=(status=true, freq=1.0))
```
"""
function ic_collapse(nel, ν, E, ρ0, l0; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(nel))d collapse problem"
    # Geometry
    dim = length(nel)
    L = dim == 2 ? [10.0, 1.25*l0] : [10.0, 10.0, 1.25*l0]
    # Simulation instructions
    instr = kwargser(:instr, kwargs; dim=dim)
    paths = set_paths(fid, info.sys.out; interactive=false)
    T0    = instr[:dtype].T0  
    T1,T2 = first(T0),last(T0) 
    L,nel = T2.(L),T1.(nel) 
    # mesh & mpts initial conditions
    ni    = T1(2)
    mesh  = setup_mesh(nel, L, instr)
    cmpr  = setup_cmpr(mesh,instr; E=T2(E), ν=T2(ν), ρ0=T2(ρ0))
    mpts    = setup_mpts(mesh, cmpr; define=geom_collapse(mesh, cmpr, ni; ℓ₀=l0))
    # time parameters
    tg    = ceil((1.0/cmpr.c)*(2.0*l0)*40.0)
    te    = 1.25*tg
    time  = setup_time(T2; te=te,tg=tg) 
    # display summary
    @info ic_log(mesh,mpts,time)
    return (;mesh, mpts, cmpr, time), (;instr, paths)
end

"""
    collapse(ic::NamedTuple, cfg::NamedTuple) -> NamedTuple

Runs the explicit solution workflow for the collapse problem using a deep copy of the initial conditions and configuration.
This function is suitable for workflows where you do not want to mutate the input data.

# Arguments
- `ic::NamedTuple`: Initial mesh/material/constitutive/time condition.
- `cfg::NamedTuple`: Simulation instructions and output paths.

# Returns
- `NamedTuple`: Simulation output with all fields from `elastoplasm` and an added `success=true` field.

# Example
```julia
result = collapse(ic, cfg)
if result.success
    println("Simulation completed successfully!")
end
```
"""
function collapse(ic::NamedTuple, cfg::NamedTuple)
    @info "Execution of collapse()"; config_plot()
    # forward-euler explicit workflow
    out = elastoplasm(deepcopy(ic), deepcopy(cfg);)
    # return output with success flag
    return out = (; out..., success=true,)
end

"""
    collapse!(ic::NamedTuple, cfg::NamedTuple) -> NamedTuple

Run the explicit solution workflow for the collapse problem, mutating the input initial conditions and configuration.
Use this when you want changes to `ic` and `cfg` to persist after the simulation.

# Arguments
- `ic::NamedTuple`: Initial mesh/material/constitutive/time condition.
- `cfg::NamedTuple`: Simulation instructions and output paths.

# Returns
- `NamedTuple`: Simulation output with all fields from `elastoplasm` and an added `success=true` field.

# Example
```julia
result = collapse!(ic, cfg)
if result.success
    println("Simulation completed successfully!")
end
```
"""
function collapse!(ic::NamedTuple, cfg::NamedTuple)
    @info "Explicit solution to collapse problem"; config_plot()
    # forward-euler explicit workflow
    out = elastoplasm(ic, cfg;)
    # return output with success flag
    return out = (; out..., success=true,)
end

#=
    plot = (;status=true,freq=1.0,what=["P"],dims=(500.0,250.0),)
    nel  = [5,10]
    # initial parameters 
    ν,E,ρ0,l0 = 0.0,1.0e4,80.0,10.0
    ic, cfg = ic_collapse(nel, ν, E, ρ0, l0; plot);
    out = collapse(ic, cfg);
=#