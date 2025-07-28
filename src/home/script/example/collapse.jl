export collapse, ic_collapse

"""
    ic_collapse(dim::Int, nel::Vector{Int64}, ν, E, ρ0, l0; fid::String=..., kwargs...) -> NamedTuple, NamedTuple

Initializes the mesh, material points, constitutive model, and simulation configuration for a collapse test.

# Arguments
- `dim::Int`: Number of spatial dimensions (2 or 3).
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `ν, E, ρ0, l0`: Physical parameters.
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration.

# Returns
- `(ic, cfg)`: Two named tuples containing mesh/material/compression info (`ic`) and instructions/paths (`cfg`).
"""
function ic_collapse(nel::Vector{Int64}, ν, E, ρ0, l0; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(nel))d collapse problem"
    # Geometry
    dim = length(nel)
    L = dim == 2 ? [10.0, 1.25*l0] : [10.0, 10.0, 1.25*l0]
    # Simulation instructions
    instr = kwargser(:instr, kwargs; dim=dim)
    paths = set_paths(fid, info.sys.out; interactive=false)
    # mesh & mp initial conditions
    ni    =  2
    mesh  = setup_mesh(nel, L, instr)
    cmpr  = setup_cmpr(dim, instr; E=E, ν=ν, ρ0=ρ0)
    mp    = setup_mps(mesh, cmpr; define=inicollapse(mesh, cmpr, ni; ℓ₀=l0))
    # time parameters
    tg    = ceil((1.0/cmpr.c)*(2.0*l0)*40.0)
    te    = 1.25*tg
    time  = (; t = [0.0, te], te = te, tg = if tg > te te else tg end, tep = nothing,)
    # display summary
    @info ic_log(mesh,mp,time)
    return (;mesh, mp, cmpr, time), (;instr, paths)
end

"""
    collapse(ic::NamedTuple, cfg::NamedTuple) -> Bool

Runs the explicit solution workflow for the collapse problem, including simulation and postprocessing.

# Arguments
- `ic::NamedTuple`: Initial mesh/material/compression/time/gravity configuration.
- `cfg::NamedTuple`: Simulation instructions and output paths.

# Returns
- `Bool`: Returns `true` if the simulation and postprocessing complete successfully.
"""
function collapse(ic::NamedTuple, cfg::NamedTuple)
    @info "Execution of collapse()"
    # extract mesh, mp, cmpr, instr, paths
    mesh,mp,cmpr = deepcopy(ic[:mesh]  ), deepcopy(ic[:mp]    ), deepcopy(ic[:cmpr])
    time         = deepcopy(ic[:time]  )
    instr,paths  = deepcopy(cfg[:instr]), deepcopy(cfg[:paths])                                              
    # action
    elastoplasm(mp,mesh,cmpr,time,paths,instr)
    # return success
    return (;sucess=true,mesh,mp)
end

#=
    plot = (;status=true,freq=1.0,what=["P"],dims=(500.0,250.0),)
    nel  = [5,10]
    # initial parameters 
    ν,E,ρ0,l0 = 0.0,1.0e4,80.0,10.0
    ic, cfg = ic_collapse(nel, ν, E, ρ0, l0; plot);
    out = collapse(ic, cfg);
=#