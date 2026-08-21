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
    L = dim == 2 ? [l0, 1.25*l0] : [2.0, 2.0, 1.25*l0]

    ly = 1.25*l0

    lx = ly/nel[2]

    L = [lx, ly]


    # init & kwargs
    instr = kwargser(kwargs; dim=dim)
    paths = set_paths(fid, self.sys.out; interactive=false)
    T0    = instr.dtype.T0  
    T1,T2 = first(T0),last(T0) 
    L,nel = T2.(L),T1.(nel) 
    # mesh & mpts initial conditions
    # mesh & mpts initial conditions
    ni    = T1(2)
    geom  = setup_geometry(L,nel,instr)
    # mesh, mpts, cmpr & time initial conditions
    mesh  = setup_mesh(geom, instr)
    cmpr  = setup_cmpr(mesh      ; E=T2(E), ν=T2(ν), ρ0=T2(ρ0))
    mpts  = setup_mpts(mesh, instr, cmpr; geom = get_collapse(mesh, cmpr, ni; ℓ₀=l0))
    # time parameters
    tg    = ceil((1.0/cmpr.c)*(2.0*l0)*40.0)
    te    = tg
    time  = setup_time(instr; te=te,tg=tg) 
    # display summary
    @info ic_log(mesh,mpts,time,instr)
    misc = (;
        prefix = "$((length(mesh.prprt.nel)-1))d_$(instr.transfer.trsfr)"
    )
    # export to jld2 file and return path
    return export_setup(mesh,mpts,cmpr,time,instr,paths,misc; path = paths[:dat], file = "collapse_simulation")
end