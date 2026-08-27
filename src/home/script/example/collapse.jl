export collapse_problem

"""
    collapse_problem(nel::Vector{Int64}; fid::String=..., kwargs...) -> String

Initializes the mesh, material points, and simulation configuration for a dry
granular Drucker-Prager block-collapse test — a plain rectangular block released
under gravity with no pre-existing slope — and exports it to a `.jld2` file (same
pipeline as `slump_problem`). Modeled directly on MaterialPointSolver.jl's
`2d_collapse.jl` reference scenario (`LandslideSIM/Archive_MaterialPointSolver.jl_paper`):
domain (`L=[0.92,0.22]`), block size (`w0=0.2,h0=0.1`), and material constants
(`ρ0=2650`, `E` derived from `Ks=7e5`, `ν=0.3`, `ϕ=19.8°`, cohesionless `c0=0`) are all
hardcoded to match that scenario directly — `nel` is the only argument that normally
needs to change (mesh resolution), everything else is a keyword override for the rare
case it's actually needed. Contrast with `slump_problem` (pre-carved sloped geometry)
and `column_problem` (1-D elastic self-weight column, no plasticity — renamed from the
old `collapse_problem`, see `src/home/script/example/column.jl`).

# Arguments
- `nel::Vector{Int64}`: Number of elements in each dimension — the one argument meant
  to vary from call to call (mesh resolution).
- `fid::String`: (Optional) File or run identifier.
- `kwargs...`: Additional keyword arguments for simulation configuration (see
  `get_solver`), plus overrides for `L`,`w0`,`h0`,`ρ0`,`E`,`ν`,`ϕ`,`c0` if ever needed —
  all hardcoded to the MaterialPointSolver.jl-matching values above by default.
  `plast.status=true`/`plast.constitutive="DP"` are set by default so the example
  demonstrates plasticity out of the box; base boundary is `:fixed` (sticky/no-slip
  floor) rather than the package default `:roller` (frictionless normal-only slip) —
  a granular collapse needs a rough/rigid floor for the material to actually pile up
  instead of sliding indefinitely, matching MaterialPointSolver.jl's own boundary
  treatment (fixes both `vx`/`vy` at the base, `vx` only at the side walls, which stay
  `:roller` here too).

# Returns
- `String`: Path to the exported `.jld2` setup file.

# Example
```julia
jld2 = collapse_problem([92, 22])
out  = elastoplasm!(jld2; workflows=[elastodynamic!, elastoplastic!])
```
"""
function collapse_problem(nel;
        L=[0.92,0.22], w0=0.2, h0=0.1,
        ρ0=2650.0, E=7.0e5*3.0*(1.0-2.0*0.3), ν=0.3, ϕ=deg2rad(19.8), c0=0.0,
        fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d granular collapse problem"
    # get solver and paths
    solver = get_solver(; dim=length(L),
        (; plast=(;status=true,constitutive="DP"),
           bcs=(;dirichlet=[:roller :roller; :fixed :roller]))..., kwargs...)
    paths  = mkpaths(fid,self.sys.dump;interactive=false)
    # mesh, mpts, mat & time initial conditions
    ni      = 2
    geom    = setup_geometry(L,nel,solver)
    mesh    = setup_mesh(geom, solver)
    mat     = setup_material_constants(solver; E=E, ν=ν, ρ0=ρ0)
    mpts    = setup_mpts(mesh, solver, mat; geom = get_collapse(mesh, mat, solver; ni=ni, w0=w0, h0=h0, ϕ=ϕ, c0=c0))
    time    = setup_time(solver; te=10.0, tg=10.0, tep=5.0)
    problem = MechanicalProblem(mesh,mpts,time)
    basis   = setup_basis(problem,solver)
    # export to jld2 file and return path
    return export_problem(problem,basis,solver,paths; file = "collapse_simulation")
end
