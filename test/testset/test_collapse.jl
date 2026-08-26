@testset "+ $(basename(@__FILE__))" verbose = true begin

# Configuration Generators
# ==============================================================================

"""
    generate_nels_cases() -> Vector{Vector{Int}}

Generate all combinations of element numbers.
"""
function generate_nels_cases()
    return [[5, 5], [5, 10], [5, 20], [5, 40],]
    #return [[5, 10], [5, 20],]
end

# Helper Functions
# ==============================================================================

"""
    load_simulation_setup(sim_file::String) -> NamedTuple

Load mesh, material points, and configuration from JLD2 simulation file.
"""
function load_simulation_setup(sim_file::String)
    jldopen(sim_file, "r") do file
        problem = file["ic"]["problem"]
        return (
            mesh   = problem.mesh,
            mpts   = problem.mpts,
            time   = problem.time,
            solver = file["cfg"]["solver"],
            paths  = file["cfg"]["paths"],
            misc   = file["cfg"]["misc"],
        )
    end
end

"""
    compute_collapse_error(sim_file, l0) -> Float64

Compute error between numeric and analytic solution for elastic collapse.
"""
function compute_collapse_error(jld2, l0, solver)
    # Load setup to get initial material point positions
    setup = load_simulation_setup(jld2)
    idx   = length(setup.mesh.prprt.nel)-1
    z0    = getindex.(setup.mpts.x, idx)

    # Run simulation
    out = elastoplasm!(jld2; workflows=[solver])

    # Reload to get final state
    jldopen(jld2, "r") do file
        problem = file["ic/problem"]
        mesh = problem.mesh
        mpts = problem.mpts
        ρ0   = mpts.s.ρ₀[1]

        # Numeric and analytic solution
        # mpts.s.σᵢⱼ holds typed CauchyStress objects (see concrete/tensor.jl), so pull the
        # Voigt view first rather than indexing the stored (p,dev) split directly
        xnum = abs.(getindex.(get_voigt.(mpts.s.σᵢⱼ), idx))
        ynum = z0
        x    = abs.(ρ0 * 9.81 * (l0 .- z0))
        y    = z0
        err  = sum(sqrt.((xnum .- x).^2) .* mpts.Ω₀) / (9.81 * ρ0 * l0 * sum(mpts.Ω₀))

        return err
    end
end

"""
    run_collapse_convergence_tests(basis, l0, plot_cfg, fwrk_cfg) -> Nothing

Run convergence tests for a given basis configuration.
"""
function run_collapse_convergence_tests(solver, fwrk)
    l0     = 50.0
    nels   = generate_nels_cases()
    nk     = length(nels) + 1
    errors = zeros(nk)
    hs     = zeros(nk)
    errors[1], hs[1] = Inf, Inf
    plot_path = ""
    for (k, nel) ∈ enumerate(nels)
        @testset "- nel = $nel" verbose = true begin
            strain   = (; deform = fwrk.deform)
            transfer = (; trsfr = fwrk.trsfr, C_pf = fwrk.C_pf, musl = fwrk.musl)
            stab     = (; locking = fwrk.locking, damping = fwrk.damping)
            sim = collapse_problem(nel, 0.0, 1.0e4, 80.0, l0; fid = "test/collapse", strain=strain, transfer=transfer, stab=stab)

            setup = load_simulation_setup(sim)
            err = compute_collapse_error(sim, l0, solver)
            errors[k+1] = err
            hs[k+1] = setup.mesh.prprt.h[end]
            plot_path = setup.paths[:plot]
            @test errors[k+1] < errors[k]
        end
    end
    return errors, hs, plot_path
end

# Main Test Execution
# ==============================================================================

all_errors = []
all_hs = []
all_labels = []
plot_path = ""



fwrks = [
    #=
    (;
        deform = "infinitesimal",
        trsfr = "std",
        C_pf = 1.0, 
        musl = false,
        locking = false,
        damping = 0.0
    ),
    =#
    (;
        deform = "finite",
        trsfr = "std",
        C_pf = 1.0, 
        musl = false,
        locking = false,
        damping = 0.0
    ),
]

solver = elastoquasistatic!

for fwrk in fwrks
    @info "Running convergence tests for workflow: $solver"
    @testset "- 2d elastic collapse" verbose = true begin
        errors, hs, sim_plot_path = run_collapse_convergence_tests(solver, fwrk)
        plot_path = sim_plot_path
        push!(all_errors, errors)
        push!(all_hs, hs)
        push!(all_labels, "$(string(solver)), $(fwrk.deform))")
    end
end

# Create convergence plot
config_plot()
opts = (;
    dims    = (450, 250),
    backend = gr(legend=true, markersize=2.0, markershape=:circle, markerstrokewidth=0.75),
    tit     = "Convergence for elastic collapse",
    file    = joinpath(plot_path, "2d_elastic_column_collapse_convergence.png"),
)
opts.backend




for j in 1:length(all_errors)
    if j == 1
        p = plot(
            1.0 ./ all_hs[j][2:end], all_errors[j][2:end],
            seriestype = :line,
            xlabel     = L"h^{-1} \,\, [\mathrm{m}^{-1}]",
            ylabel     = L"\epsilon" * " [-]",
            yscale     = :log10,
            label      = all_labels[j],
            title      = opts.tit,
            size       = opts.dims,
        )
    else
        p = plot!(
            1.0 ./ all_hs[j][2:end], all_errors[j][2:end],
            seriestype = :line,
            label      = all_labels[j],
            title      = opts.tit,
            size       = opts.dims,
        )
    end
end
display(current())
save_plot(opts)

end # @testset