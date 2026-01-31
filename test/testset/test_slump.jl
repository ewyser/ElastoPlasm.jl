# Test configuration constants
const DEFORMATIONS = ["finite", "infinitesimal"]
const TRANSFERS = ["std", "tpic", "apic"]
const LOCKING_OPTIONS = [true, false]

const TEST_GEOMETRIES = [
    (dim=2, L=[64.1584, 64.1584/4.0], nel=[40, 10]),
    (dim=3, L=[64.1584, 64.1584/4.0, 64.1584/4.0], nel=[40, 10, 10]),
]

const BASIS_CONFIGS = [
    (which="bsmpm", how=nothing, ghost=false),
    #(which="gimpm", how="undeformed", ghost=true),
    #(which="smpm", how=nothing, ghost=true),
]

"""
    generate_test_cases() -> Vector{NamedTuple}

Generate all combinations of deformation frameworks, transfer schemes, and locking options.
"""
function generate_test_cases()
    return [(deform=d, trsfr=t, locking=l) 
            for d in DEFORMATIONS 
            for t in TRANSFERS 
            for l in LOCKING_OPTIONS]
end

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
    get_test_plot_config() -> NamedTuple

Generate plot configuration for slump tests using new plot format.
"""
function get_test_plot_config()
    return (
        status = true,
        dpi    = 500,
        freq   = 1.0,
        what   = [(;mpts=(name="epII", cblim=(0.0, 1.5)))],
        dims   = (500.0, 250.0),
    )
end

"""
    get_test_grf_config() -> NamedTuple

Generate Gaussian Random Field configuration for tests.
"""
function get_test_grf_config()
    return (
        status      = true,
        covariance  = "gaussian",
        param       = (Iₓ=[2.5, 2.5, 2.5], Nₕ=5000, kₘ=100),
    )
end

"""
    run_slump_test_case(setup, test_case, basis, plot_cfg, grf_cfg, case_id) -> Bool

Execute a single slump test configuration and return success status.
"""
function run_slump_test_case(setup, test_case, basis, plot_cfg, grf_cfg, case_id)
    kwargs = Dict(
        :basis   => basis,
        :fwrk    => (
            deform  = test_case.deform,
            trsfr   = test_case.trsfr,
            musl    = true,
            C_pf    = 0.99,
            locking = test_case.locking,
            damping = 0.1,
        ),
        :nonloc  => (status=false, ls=0.5),
        :plot    => plot_cfg,
        :grf     => grf_cfg,
        :backend => setup.instr.backend,
    )
    
    instr = kwargser(kwargs; dim=setup.mesh.prprt.dim)
    misc = (prefix = "case_$(case_id)_$(setup.mesh.prprt.dim)d_$(instr.fwrk.trsfr)_$(instr.fwrk.deform)",)
    
    sim = export_setup(
        setup.mesh, setup.mpts, setup.cmpr, setup.time, 
        instr, setup.paths, misc;
        path = setup.paths[:dat],
        file = "slump_simulation"
    )
    
    status = false
    try
        @suppress begin
            status = elastoplasm(sim; workflow=[elastodynamic!, elastoplastic!]).success
        end
    catch e
        @warn "Test case failed" test_case exception=(e, catch_backtrace())
    end
    
    return status
end

"""
    run_slump_test_suite(sim_file, basis, geometries; progress_desc="Testing") -> Nothing

Run complete test suite for given simulation file and basis configuration.
"""
function run_slump_test_suite(sim_file, basis, geometry; progress_desc="Testing")
    setup = load_simulation_setup(sim_file)
    test_cases = generate_test_cases()
    plot_cfg = get_test_plot_config()
    grf_cfg = get_test_grf_config()
    
    prog = Progress(length(test_cases); dt=0.1, desc=progress_desc, barlen=10)
    
    for (k, test_case) in enumerate(test_cases)
        @testset "$(test_case.deform) + $(test_case.trsfr) + locking=$(test_case.locking)" begin
            status = run_slump_test_case(setup, test_case, basis, plot_cfg, grf_cfg, k)
            @test status
        end
        next!(prog; desc="$progress_desc $(k)/$(length(test_cases))...")
    end
    
    finish!(prog)
end

@testset "+ $(basename(@__FILE__))" verbose = true begin
    plot_cfg = get_test_plot_config()
    grf_cfg = get_test_grf_config()
    
    for basis in BASIS_CONFIGS
        @info "Testing with $(basis.which) basis"
        
        for geom in TEST_GEOMETRIES
            sim = ic_slump(geom.L, geom.nel; fid="test/slump", grf=grf_cfg)
            
            @testset "- $(geom.dim)d geometry with $(basis.which) basis" verbose = true begin
                run_slump_test_suite(
                    sim, basis, geom;
                    progress_desc="$(geom.dim)d completion:"
                )
            end
        end
    end
end