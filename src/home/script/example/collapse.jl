export collapse_problem

"""
    collapse_problem(nel::Vector{Int64}, ν, E, ρ0, l0; fid::String=..., kwargs...) -> String

Initializes the mesh, material points, and simulation configuration for a column
collapse test, and exports it to a `.jld2` file (same pipeline as `slump_problem`).

# Arguments
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `ν`: Poisson's ratio.
- `E`: Young's modulus.
- `ρ0`: Initial density.
- `l0`: Characteristic length.
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration.

# Returns
- `String`: Path to the exported `.jld2` setup file.

# Example
```julia
jld2 = collapse_problem([5, 10], 0.0, 1.0e4, 80.0, 10.0; plot = (; status=true, freq=1.0, what=["P"], dims=(500.0,250.0) ));
```
"""
function collapse_problem(nel, ν, E, ρ0, l0; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(nel))d collapse problem"
    # Geometry
    dim = length(nel)
    ly  = 1.25*l0
    lx  = ly/nel[2]
    L   = [lx, ly]
    # get solver and paths
    solver = get_solver(; dim=dim, kwargs...)
    paths  = mkpaths(fid, self.sys.dump; interactive=false)
    # mesh, mpts, mat & time initial conditions
    ni      = 2
    geom    = setup_geometry(L,nel,solver)
    mesh    = setup_mesh(geom, solver)
    mat     = setup_material_constants(solver; E=E, ν=ν, ρ0=ρ0)
    mpts    = setup_mpts(mesh, solver, mat; geom = get_collapse(mesh, mat, ni; ℓ₀=l0))
    # time parameters
    tg      = ceil((1.0/mat.c)*(2.0*l0)*40.0)
    time    = setup_time(solver; te=tg, tg=tg)
    problem = MechanicalProblem(mesh,mpts,time)
    basis   = setup_basis(problem,solver)
    # export to jld2 file and return path
    return export_problem(problem,basis,solver,paths; file = "collapse_simulation")
end
