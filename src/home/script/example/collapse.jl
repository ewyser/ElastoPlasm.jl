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
ic, cfg = ic_collapse([5, 10], 0.0, 1.0e4, 80.0, 10.0; plot = (; status=true, freq=1.0, what=["P"], dims=(500.0,250.0) ));
```
"""
function ic_collapse(nel, ν, E, ρ0, l0; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(nel))d collapse problem"
    # Geometry
    dim = length(nel)
    L = dim == 2 ? [10.0, 1.25*l0] : [10.0, 10.0, 1.25*l0]
    # Simulation instructions
    instr = kwargser(:instr, kwargs; dim=dim)
    instr = merge(instr, (bcs = (;dirichlet=[:roller :roller;:fixed :roller]),))

    paths = set_paths(fid, info.sys.out; interactive=false)
    T0    = instr[:dtype].T0  
    T1,T2 = first(T0),last(T0) 
    L,nel = T2.(L),T1.(nel) 
    # mesh & mpts initial conditions
    h       = [L[end]/nel[end],L[end]/nel[end]]
    nel     = [5,nel[end]]
    L       = [5*h[end],L[end]]
    ndim,nn = length(L),4^length(L)
    ni      = T1(2)

    geom  = (; ndim = T1(ndim), nn = T1(nn), h =T2.(h), nel = T1.(nel), L = T2.(L))
    mesh  = setup_mesh(instr     ; geom = geom     )
    cmpr  = setup_cmpr(mesh      ; E=T2(E), ν=T2(ν), ρ0=T2(ρ0))
    mpts  = setup_mpts(mesh, cmpr; geom = get_collapse(mesh, cmpr, ni; ℓ₀=l0))
    # time parameters
    tg    = ceil((1.0/cmpr.c)*(2.0*l0)*40.0)
    te    = tg
    time  = setup_time(instr; te=te,tg=tg) 
    # display summary
    @info ic_log(mesh,mpts,time,instr)
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
result = collapse(ic, cfg);
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