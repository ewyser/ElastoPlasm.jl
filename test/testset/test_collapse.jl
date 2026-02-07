@testset "+ $(basename(@__FILE__))" verbose = true begin

# Configuration Generators
# ==============================================================================

"""
    generate_basis_cases() -> Vector{NamedTuple}

Generate all combinations of shape functions.
"""
function generate_basis_cases()
    return [(which="bsmpm", how=nothing, ghost=false),
    (; which = "gimpm", how = "undeformed", ghost = true ),
    #(which="smpm", how=nothing, ghost=true),
    ]
end

"""
    generate_nels_cases() -> Vector{Vector{Int}}

Generate all combinations of element numbers.
"""
function generate_nels_cases()
    return [[5, 5], [5, 10], [5, 20], [5, 40], [5, 80]]
end

"""
    generate_fwrk_cases() -> Vector{NamedTuple}

Generate all combinations of deformation frameworks, transfer schemes, and locking options.
"""
function generate_fwrk_cases()
    return [(deform=d, trsfr=t, locking=l, C_pf = 1.0, musl = m, damping = 0.1) 
            for d in ["finite", "infinitesimal"] 
            for t in ["std", "tpic", "apic"] 
            for l in [true, false]
            for m in [true, false]
        ]
end

# Helper Functions
# ==============================================================================

"""
    load_simulation_setup(sim_file::String) -> NamedTuple

Load mesh, material points, and configuration from JLD2 simulation file.
"""
function load_simulation_setup(sim_file::String)
    jldopen(sim_file, "r") do file
        return (
            mesh  = file["ic"]["mesh"],
            mpts  = file["ic"]["mpts"],
            cmpr  = file["ic"]["cmpr"],
            time  = file["ic"]["time"],
            instr = file["cfg"]["instr"],
            paths = file["cfg"]["paths"],
            misc  = file["cfg"]["misc"],
        )
    end
end

"""
    compute_collapse_error(sim_file, l0) -> Float64

Compute error between numeric and analytic solution for elastic collapse.
"""
function compute_collapse_error(jld2, l0)
    # Load setup to get initial material point positions
    setup = load_simulation_setup(jld2)
    z0 = copy(setup.mpts.x[end, :])
    
    # Run simulation
    out = elastoplasm!(jld2; workflow=[elastodynamic!])
    
    # Reload to get final state
    jldopen(jld2, "r") do file
        mesh = file["ic/mesh"]
        mpts = file["ic/mpts"]
        cmpr = file["ic/cmpr"]
        
        # Numeric and analytic solution
        idx  = mesh.prprt.dim == 2 ? 2 : 3
        xnum = abs.(mpts.s.σᵢ[idx, :])
        ynum = z0
        x    = abs.(cmpr.ρ0 * 9.81 * (l0 .- z0))
        y    = z0
        err  = sum(sqrt.((xnum .- x).^2) .* mpts.Ω₀) / (9.81 * cmpr.ρ0 * l0 * sum(mpts.Ω₀))
        
        return err
    end
end

"""
    run_collapse_convergence_tests(basis, l0, plot_cfg, fwrk_cfg) -> Nothing

Run convergence tests for a given basis configuration.
"""
function run_collapse_convergence_tests(basis, l0, plot_cfg, fwrk_cfg)
    nels   = generate_nels_cases()
    nk     = length(nels) + 1
    errors = zeros(nk)
    hs     = zeros(nk)
    errors[1], hs[1] = Inf, Inf
    
    for (k, nel) ∈ enumerate(nels)
        @testset "- nel = $nel" verbose = true begin
            sim = ic_collapse(nel, 0.0, 1.0e4, 80.0, l0; 
                             fid = "test/collapse", 
                             plot = plot_cfg, 
                             basis = basis, 
                             fwrk = fwrk_cfg)
            
            setup = load_simulation_setup(sim)
            err = compute_collapse_error(sim, l0)
            errors[k+1] = err
            hs[k+1] = setup.mesh.prprt.h[end]
            @test errors[k+1] < errors[k]
        end
    end
    
    # Return errors and mesh sizes for convergence plotting
    return errors, hs
end

# Main Test Execution
# ==============================================================================

l0 = 50.0
plot_cfg = get_test_plot_config()
fwrk_cfg = get_test_fwrk_config()

# Storage for convergence plotting
all_errors = []
all_hs = []
all_labels = []
plot_path = ""

for (j, basis) in enumerate(BASIS_CONFIGS)
    @testset "- 2d geometry with $(basis.which) basis" verbose = true begin
        errors, hs = run_collapse_convergence_tests(basis, l0, plot_cfg, fwrk_cfg)
        push!(all_errors, errors)
        push!(all_hs, hs)
        push!(all_labels, basis.which)
        
        # Get plot path from first test case
        if j == 1 && length(TEST_NELS) > 0
            sim = ic_collapse(TEST_NELS[1], 0.0, 1.0e4, 80.0, l0; 
                             fid = "test/collapse", 
                             plot = plot_cfg, 
                             basis = basis, 
                             fwrk = fwrk_cfg)
            setup = load_simulation_setup(sim)
            plot_path = setup.paths[:plot]
        end
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

for (j, basis) in enumerate(BASIS_CONFIGS)
    get_plot = j == 1 ? plot : plot!
    p = get_plot(
        1.0 ./ all_hs[j][2:end], all_errors[j][2:end],
        seriestype = :line,
        xlabel     = L"h^{-1} \,\, [\mathrm{m}^{-1}]",
        ylabel     = L"\epsilon" * " [-]",
        yscale     = :log10,
        label      = "2D $(all_labels[j]), $(fwrk_cfg.trsfr) mapping",
        title      = opts.tit,
        size       = opts.dims,
    )
end
display(current())
save_plot(opts)

end # @testset