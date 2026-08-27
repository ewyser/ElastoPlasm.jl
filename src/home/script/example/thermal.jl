export thermal_problem

"""
    thermal_problem(L::Vector{Float64}, nel::Vector{Int64}; fid::String=..., kwargs...) -> String

Initializes the mesh, material points, and simulation configuration for a heat-conduction
test problem, and exports it to a `.jld2` file (same pipeline as `slump_problem`).

# Arguments
- `L::Vector{Float64}`: Domain dimensions.
- `nel::Vector{Int64}`: Number of elements in each dimension.
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration (see `get_solver`);
  override the thermal-specific defaults below the same way `slump_problem` does
  (`thermal_problem(L, nel; cli()..., basis=(;get_default().basis...,which="gimpm",how="Uii"))`
  — spread `get_default().basis` first since `basis` now has more fields than just
  `which`/`how`, per `get_solver`'s shallow per-key merge).

# Returns
- `String`: Path to the exported `.jld2` setup file.

# Example
```julia
jld2 = thermal_problem([16.0, 16.0], [80, 80])
out  = elastoplasm(jld2; workflows=[thermodynamic!])
```
"""
function thermal_problem(L,nel; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d thermal problem"
    # get solver and paths
    solver = get_solver(; dim=length(L),
        (;
            basis = (; get_default().basis..., which = "bsmpm", how = nothing),
            bcs   = (; dirichlet = fill(:fixed, length(L), 2)), # one [lower upper] row per dimension
            plot  = (; status = true, freq = 1.0, dpi = 300, what = [(;mpts=get_mpts_variable_config()["T"])]),
        )...,
        kwargs...)
    paths   = mkpaths(fid,self.sys.dump;interactive=false)
    # mesh, mpts, mat & time initial conditions
    geom    = setup_geometry(L,nel,solver)
    problem = setup_problem(geom,solver,get_thermal; te = 5.0, tg = 5.0, tep = 0.0, thermal = true)
    basis   = setup_basis(problem,solver)
    # export to jld2 file and return path
    return export_problem(problem,basis,solver,paths; file = "thermal_simulation")
end
