@testset "+ $(basename(@__FILE__))" verbose = true begin
    EP = ElastoPlasm
    @testset "shapes" begin
        @test EP.inside(Block([0, 0], [1, 1]), [1.0, 0.0]) && !EP.inside(Block([0, 0], [1, 1]), [1.0 + 1e-12, 0.5])   # bounds included
        @test EP.inside(Block([0, 0], [1, 1]), [0.5, 9.0, 1.0])                                    # (x, z) box extruded along y
        @test EP.inside(Block([0, 0, 0], [1, 1, 1]), [0.5, 0.5, 0.5]) && !EP.inside(Block([0, 0, 0], [1, 1, 1]), [0.5, 2.0, 0.5])
        @test EP.inside(Ball([0, 0], 1.0), [0.6, 0.8]) && !EP.inside(Ball([0, 0], 1.0), [0.6, 0.81])
        @test EP.inside(Ball([0, 0, 0], 1.0), [0.0, 0.6, 0.8]) && EP.inside(Ball([0, 0], 1.0), [0.6, 5.0, 0.8])   # sphere; cylinder
        ell = Polygon([[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]])                          # concave
        @test EP.inside(ell, [0.5, 1.5]) && EP.inside(ell, [1.5, 0.5]) && !EP.inside(ell, [1.5, 1.5])
        @test EP.inside(Region(x -> x[1] > x[end]), [2.0, 1.0])
        @test_throws ArgumentError Polygon([[0, 0], [1, 0]])
        @test_throws ArgumentError Block([1, 0], [0, 1])
    end
    @testset "specification and validation" begin
        s = Case(; L=[4, 2], nel=[4, 2], zones=[Zone("a", Block([0, 0], [4, 2]))], material=Dict("c0" => 1.0e3),
                   solver=(; basis=(; which="gimpm"), bcs=(; dirichlet=[:fixed :fixed; :fixed :free])))
        @test s.material.c0 == 1.0e3 && s.material.E == MATERIAL.E && s.bcs == [:fixed :fixed; :fixed :free] && !haskey(s.solver, :bcs)
        @test s.run == EP.RUN && EP.case_solver_options(s).solution == "explicit"
        q = Case(; L=[4, 2], nel=[4, 2], width=1.0, nely=2, run=(; workflow="elastoquasistatic!"))
        @test EP.domain(q) == ([4.0, 1.0, 2.0], [4, 2, 2]) && EP.case_solver_options(q).solution == "implicit"
        @test EP.case_solver_options(q).bcs.dirichlet == fill(:roller, 3, 2)
        @test_throws ArgumentError Case(; L=[1, 1], nel=[1, 1], material=(; Ee=1.0))
        @test_throws ArgumentError Case(; L=[1, 1], nel=[1, 1], solver=(; basis=(; whch="x")))
        @test_throws ArgumentError Case(; L=[1, 1], nel=[1, 1], run=(; workflow="x"))
        @test_throws ArgumentError Case(; L=[1, 1], nel=[1, 1], bcs=[:fixed :glued; :free :free])
        @test_throws ArgumentError Case(; L=[1, 1, 1], nel=[1, 1, 1], width=1.0, nely=1)
        @test_throws ArgumentError Zone("z", Block([0, 0], [1, 1]); EE=1.0)
    end
    @testset "points: last zone wins, overrides, velocity" begin
        s = Case(; L=[4, 2], nel=[4, 2], ni=2,
                   zones=[Zone("all", Block([0, 0], [4, 2])), Zone("left", Block([0, 0], [1, 2]); E=5e6, phi=0.1, v=[1.0, -2.0]),
                          Zone("none", Block([10, 10], [11, 11]))])
        solver = get_solver(; dim=2, plot=(; status=false))
        mesh   = EP.setup_mesh(s.L, s.nel, solver)
        p      = case_points(s, mesh, solver)
        left   = p.xp[1, :] .< 1
        @test size(p.xp, 2) == 32 && count(left) == 8 && all(p.zone[left] .== 2) && all(p.zone[.!left] .== 1)
        @test all(p.E[left] .== 5e6) && all(p.E[.!left] .== MATERIAL.E) && all(p.phi[left] .== 0.1) && all(p.c0 .== MATERIAL.c0)
        @test all(p.vp[:, left] .== [1.0, -2.0]) && all(p.vp[:, .!left] .== 0)
        x3 = case_points(Case(; L=[4, 2], nel=[4, 2], width=1.0, nely=1, zones=s.zones), 
                         (sv = get_solver(; dim=3, plot=(; status=false)); EP.setup_mesh([4.0, 1.0, 2.0], [4, 1, 2], sv)), sv)
        @test size(x3.xp, 2) == 64 && count(x3.zone .== 2) == 16 && all(x3.vp[:, x3.zone .== 2] .== [1.0, 0.0, -2.0])   # (x, z) velocity
    end
    @testset "TOML round trip" begin
        s = Case(; L=[20, 10], nel=[20, 10], ni=3, width=4.0, nely=4, bcs=[:roller :roller; :fixed :free], thermal=true,
                   zones=[Zone("soil", Polygon([[0, 0], [20, 0], [20, 6], [0, 6]])), Zone("rock", Block([0, 0], [20, 2]); E=5e7, phi=30.0*π/180),
                          Zone("boulder", Ball([10, 4], 1.0); rho0=3000.0, v=[0.0, 0.0, -1.0])],
                   material=(; c0=1.0e4, phi=25.0*π/180), solver=(; basis=(; which="gimpm", C_pf=0.37), stab=(; locking=false),
                   plot=(; what=EP.option_value(["plot", "what"], ["pressure"]))), run=(; folder="x/y", te=1.0, workflow="elastodynamic!"))
        back = load_case(save_case(joinpath(mktempdir(), "p.toml"), s))
        @test (back.L, back.nel, back.ni, back.width, back.nely, back.bcs, back.thermal, back.run) == (s.L, s.nel, s.ni, s.width, s.nely, s.bcs, s.thermal, s.run)
        @test [z.name for z ∈ back.zones] == ["soil", "rock", "boulder"] && back.zones[3].shape.radius == 1.0 && back.zones[3].props[:v] == [0.0, 0.0, -1.0]
        @test back.zones[2].props[:phi] ≈ 30.0*π/180 rtol=1e-12
        @test back.material.phi ≈ s.material.phi rtol=1e-12
        @test back.material.c0 == 1.0e4 && Dict(pairs(back.solver.basis)) == Dict(:which => "gimpm", :C_pf => 0.37)
        @test back.solver.stab == (; locking=false) && back.solver.plot.what[1].mpts.key == "P"
        @test_throws ArgumentError save_case(joinpath(mktempdir(), "r.toml"), Case(; L=[1, 1], nel=[1, 1], zones=[Zone("f", Region(x -> true))]))
    end
    @testset "checks and time step" begin
        layered(; kw...) = Case(; L=[20, 10], nel=[20, 10], kw...,
            zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]); E=5e7, nu=0.2, rho0=2900.0)])
        s = layered(; solver=(; plot=(; status=false)), run=(; folder="test/spec_dt"))
        c = check_case(s)
        @test isempty(c.errors) && isempty(c.warnings) && startswith(c.info[1], "480 material points, 200 elements (2D)")
        pb = load_simulation(setup_simulation(s)).problem
        @test pb.mpts.nmp == 480
        @test estimate_dt(s)[1] ≈ EP.get_dt(pb.mpts, pb.mesh.prprt, pb.time, Inf) rtol=1e-12   # = the solver's
        ex = layered(; width=4.0, nely=4, solver=(; plot=(; status=false)), run=(; folder="test/spec_dt3"))
        pb3 = load_simulation(setup_simulation(ex)).problem
        @test startswith(check_case(ex).info[1], "$(pb3.mpts.nmp) material points, 800 elements (3D") && pb3.mpts.nmp == 480*8
        @test estimate_dt(ex)[1] ≈ EP.get_dt(pb3.mpts, pb3.mesh.prprt, pb3.time, Inf) rtol=1e-12
        none = Case(; L=[20, 10], nel=[20, 10], zones=[Zone("away", Polygon([[30, 0], [40, 0], [40, 5]]))])
        @test length(check_case(none).errors) == 1 && length(check_case(none).warnings) == 2   # no points; outside; empty zone
        hidden = Case(; L=[20, 10], nel=[20, 10], zones=reverse(s.zones))
        @test length(check_case(hidden).warnings) == 1                                          # rock covered by soil
        @test length(check_case(layered(; solver=(; material=(; elastic="hypoelastic")), run=(; workflow="elastoquasistatic!"))).errors) == 1
        @test length(check_case(layered(; solver=(; solution="implicit"))).errors) == 1
        @test length(check_case(layered(; solver=(; dtype=(; precision="32-bit precision (or single precision)"), basis=(; which="gimpm")))).warnings) == 1
        @test length(check_case(layered(; run=(; tg=20.0))).warnings) == 1
        @test length(check_case(Case(; L=[20, 10], nel=[10, 10], zones=s.zones, run=(; workflow="elastoquasistatic!"))).warnings) == 1
    end
    @testset "zone constants, boundaries and options reach the problem" begin
        s = Case(; L=[20, 10], nel=[20, 10], bcs=[:roller :roller; :fixed :free],
                   zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]); E=5e7, nu=0.2, rho0=2900.0)],
                   solver=(; basis=(; which="gimpm", C_pf=0.37), stab=(; locking=false), material=(; plastic="VM"), plot=(; status=false)),
                   run=(; folder="test/case_zones"))
        pb = load_simulation(setup_simulation(s)).problem; solver = pb.solver
        m = pb.mpts
        @test solver.basis == (; which="gimpm", trsfr=get_default().basis.trsfr, C_pf=0.37) && !solver.stab.locking &&
              solver.stab.damping == get_default().stab.damping && solver.material.plastic == "VM" && solver.bcs.dirichlet == s.bcs
        rock = getindex.(m.x, 2) .< 2.0
        K, G = EP.get_elastic_moduli(5e7, 0.2); K0, G0 = EP.get_elastic_moduli(MATERIAL.E, MATERIAL.nu)
        @test all(c -> (c.Kc, c.Gc) == (K, G), m.s.cmp[rock]) && all(m.s.ρ₀[rock] .== 2900.0)
        @test all(c -> (c.Kc, c.Gc) == (K0, G0), m.s.cmp[.!rock]) && all(m.s.ρ₀[.!rock] .== MATERIAL.rho0)
        base = [i for i ∈ eachindex(pb.mesh.x) if pb.mesh.x[i][2] == 0]
        @test all(pb.mesh.s.bcs.status[:, base])
    end
    @testset "free option values" begin
        @test EP.option_value(["basis", "C_pf"], 0.37) == 0.37 && EP.option_value(["grf", "param", "Nₕ"], "2000") === 2000
        @test EP.option_value(["grf", "param", "Iₓ"], "3, 4, 5") == [3.0, 4.0, 5.0] && EP.option_value(["stab", "locking"], "false") === false
        @test_throws ArgumentError EP.option_value(["grf", "param", "Nₕ"], 2.5)                    # an integer option
        @test_throws ArgumentError EP.option_value(["basis", "which"], "foo")                      # not a numeric option
        @test [v.mpts.key for v ∈ EP.option_value(["plot", "what"], ["pressure", "plastic strain"])] == ["P", "epII"]
    end
    @testset "simulation: build, save, load, execute!" begin
        dir = joinpath(mktempdir(), "mycase")
        c = Case(; L=[4, 2], nel=[4, 2], zones=[Zone("a", Block([0, 0], [4, 1]))], solver=(; plot=(; status=false, freq=0.1)),
                   run=(; folder=dir * "/", te=0.2, tg=0.1, tep=0.0, workflow="elastodynamic!"))
        sim = Simulation(Problem(c); c.run.workflow, c.run.lstps, c.run.folder)
        @test sim.folder == dir && !isdir(dir)                                              # built, nothing written
        @test Simulation(sim.problem; folder="a/b").folder == joinpath(EP.self.sys.dump, "a", "b")   # a name: under dump/
        @test_throws ArgumentError Simulation(sim.problem; workflow="nope")
        @test_throws ArgumentError Simulation(sim.problem; lstps=0)
        jld2 = setup_simulation(c)
        @test jld2 == joinpath(dir, "mycase_simulation.jld2") && isfile(joinpath(dir, "mycase_case.toml"))
        l = load_simulation(jld2); @test (l.workflow, l.lstps, l.folder) == (sim.workflow, sim.lstps, sim.folder)
        state(p) = (p.mpts.x, getfield.(p.mpts.s.σᵢ, :p), p.mpts.s.v)
        s0 = state(sim.problem); r = execute(jld2)
        @test state(r.problem) != s0 && state(load_simulation(jld2).problem) == s0               # ran a copy, file untouched
        @test execute!(jld2) == jld2 && state(load_simulation(jld2).problem) == state(r.problem)  # written back
        @test state(execute(sim).problem) == state(r.problem) && state(sim.problem) == s0          # a copy of sim
        @test state(execute!(sim).problem) == state(r.problem)                                   # in memory: the same
        seen = Bool[]; probe(p) = push!(seen, state(load_simulation(jld2).problem) == state(p))
        @test execute!(jld2; workflows=[elastodynamic!, probe]) == jld2 && seen == [true]          # other workflows, saved after each
        moved = joinpath(mktempdir(), "copied"); cp(dir, moved)
        @test load_simulation(joinpath(moved, "mycase_simulation.jld2")).folder == moved          # its folder: the file's
        @test EP.case_run((; fid="x")).folder == "x" && EP.case_run(Dict("fid" => "x", "dir" => "/a")).folder == "/a"   # older cases
        old = joinpath(mktempdir(), "old.jld2"); jldsave(old; problem=sim.problem)
        @test_throws ArgumentError load_simulation(old)                                     # an older file
        q  = Problem(EP.block_case([8, 4]; plot=(; status=false, freq=1.0, dpi=100, what=[])))   # quasi-static: implicit
        @test Simulation(q).workflow == "elastoquasistatic!"                                      # its default workflow
        @test_throws ArgumentError Simulation(q; workflow="elastodynamic!")                       # an explicit one: refused
        qs = Simulation(q; lstps=3)
        @test length(only(case_workflows(qs))(qs.problem)) == 3                                   # 3 load steps
    end
    @testset "workflows and a short run" begin
        s = Case(; L=[20, 10], nel=[20, 10], width=4.0, nely=4, zones=[Zone("soil", Block([0, 0], [20, 6]))],
                   solver=(; plot=(; status=false)), run=(; folder="test/case_run", te=0.2, tg=0.1, tep=0.0, workflow="elastodynamic!"))
        @test case_workflows(s.run.workflow, s.run.lstps) == Function[elastodynamic!]
        @test case_workflows(EP.RUN.workflow, 1) == Function[elastodynamic!, elastoplastic!]
        @test length(case_workflows("elastoquasistatic! (plastic)", 10)) == 1
        @test isfile(execute!(setup_simulation(s)))                                               # extruded, its own times
    end
end
