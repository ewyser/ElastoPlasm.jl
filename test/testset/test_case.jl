@testset "+ $(basename(@__FILE__))" verbose = true begin
    # a case written in code: the layout, then its settings in ElastoPlasm's terms, read at once (its errors here)
    function mkcase(; L, nel, ni=2, zones=Zone[], kwargs...)
        c = Case(Layout(; L, nel, ni, zones), case_settings(; kwargs...))
        case_physics(c)
        return c
    end
    EP = ElastoPlasm
    @testset "specification and validation" begin
        s = mkcase(; L=[4, 2], nel=[4, 2], zones=[Zone("a", Block([0, 0], [4, 2]))], material=Dict("c0" => 1.0e3),
                   solver=(; basis=(; which="gimpm"), bcs=(; dirichlet=[:fixed :fixed; :fixed :free])))
        @test case_physics(s).material.c0 == 1.0e3 && case_physics(s).material.E == get_material().E && case_physics(s).bcs == [:fixed :fixed; :fixed :free] && !haskey(s.settings.options, "bcs")
        @test (case_physics(s).laws, case_physics(s).grf) == (get_default().material, get_default().grf)                # the physics: defaults
        @test Solver(s).run == EP.solver_run((;)) && Solver(s).solution == "explicit"
        q = mkcase(; L=[4, 1, 2], nel=[4, 2, 2], run=(; workflow="elastoquasistatic!"))
        @test (q.layout.L, q.layout.nel) == ([4.0, 1.0, 2.0], [4, 2, 2]) && Solver(q).solution == "implicit"   # follows the workflow
        @test case_physics(q).bcs == fill(:roller, 3, 2)
        @test_throws ArgumentError mkcase(; L=[1, 1], nel=[1, 1], material=(; Ee=1.0))
        @test_throws ArgumentError mkcase(; L=[1, 1], nel=[1, 1], solver=(; basis=(; whch="x")))
        @test_throws ArgumentError mkcase(; L=[1, 1], nel=[1, 1], run=(; workflow="x"))
        @test_throws ArgumentError mkcase(; L=[1, 1], nel=[1, 1], bcs=[:fixed :glued; :free :free])
        @test_throws ArgumentError mkcase(; L=[1, 1, 1, 1], nel=[1, 1, 1, 1])                    # 2D or 3D
        z = [Zone("z", Block([0, 0], [1, 1]))]
        @test_throws ArgumentError mkcase(; L=[1, 1], nel=[1, 1], zones=z, properties=Dict("z" => (; EE=1.0)))   # get_zone_properties() only
        @test_throws ArgumentError mkcase(; L=[1, 1], nel=[1, 1], zones=z, properties=Dict("z" => (; E="stiff")))
        @test_throws ArgumentError mkcase(; L=[1, 1], nel=[1, 1], zones=z, properties=Dict("cliff" => (; E=1.0)))   # no such zone
        one = Layout(; L=[1, 1], nel=[1, 1])
        @test_throws ArgumentError case_physics(Case(one, Settings(; defaults=Dict("eta" => 1.0))))   # not a material constant
        @test_throws ArgumentError case_physics(Case(one, Settings(; sides=Dict("top" => "glued"))))
    end
    @testset "a case reads its layout and settings" begin
        s = mkcase(; L=[4, 2], nel=[4, 2], zones=[Zone("a", Block([0, 0], [4, 1]))], properties=Dict("a" => (; E=5e6)),
                   material=(; c0=1.0e3), thermal=true, bcs=[:roller :roller; :fixed :free], solver=(; basis=(; which="gimpm")),
                   run=(; te=1.0))
        l, t = s.layout, s.settings
        @test fieldnames(typeof(l.zones[1])) == (:name, :shape) && case_physics(s).properties == Dict("a" => (; E=5e6))   # geometry; values by name
        @test t.sides == Dict("bottom" => "fixed", "top" => "free")                         # a roller side: not written
        @test t.defaults == Dict("c0" => 1.0e3) && t.options == Dict("basis" => Dict("which" => "gimpm")) && t.run["thermal"]
        r = Case(l, t)                                                                       # the two alone: the same case
        @test case_physics(r) == case_physics(s)
        t.sides["bottom"] = "free"; push!(l.zones, Zone("b", Block([0, 0], [1, 1])))         # change what was given
        @test case_physics(r).bcs[2, 1] == :fixed && r.settings.sides["bottom"] == "fixed" && length(r.layout.zones) == 1   # the case: its own copy
        q = mkcase(; L=[1, 2, 3], nel=[1, 2, 3], bcs=[:free :fixed; :roller :fixed; :fixed :free])   # 3D: y is front and back
        @test q.settings.sides == Dict("left" => "free", "right" => "fixed", "back" => "fixed", "bottom" => "fixed", "top" => "free")
        @test case_physics(q).bcs == [:free :fixed; :roller :fixed; :fixed :free]
        @test EP.case_folder(Case(Layout(; L=[1, 1], nel=[1, 1], folder="/a/b"), Settings())) == "/a/b"   # a file's folder
        @test EP.case_folder(Case(Layout(; L=[1, 1], nel=[1, 1], folder="/a/b"), Settings(; run=Dict("folder" => "c")))) == "c"   # unless set
    end
    @testset "a solver solves a case" begin
        s = mkcase(; L=[8, 4], nel=[8, 4], zones=[Zone("a", Block([0, 0], [8, 2]))], properties=Dict("a" => (; E=5e6)),
                   solver=(; material=(; plastic="VM"), stab=(; damping=0.2), plot=(; status=false)), run=(; te=1.0))
        @test case_physics(s).laws == (; plastic="VM", elastic=get_default().material.elastic)
        @test fieldnames(Case) == (:layout, :settings)                                     # EasyDraw's: the values read by case_physics
        a = Solver(s)                                                                        # its numerics: from the settings
        @test a isa Solver{Int64,Float64,2} && isempty(intersect(fieldnames(typeof(a)), (:case, :material, :grf, :bcs)))
        @test a.stab.damping == 0.2 && a.run.te == 1.0
        b = Solver(s; basis=(; which="gimpm"), stab=(; locking=false), run=(; te=2.0))      # overridden, field by field
        @test (b.basis.which, b.stab.locking, b.stab.damping, b.run.te) == ("gimpm", false, 0.2, 2.0)
        @test !haskey(s.settings.options, "basis")                                           # the case: unchanged
        @test_throws ArgumentError Solver(s; material=(; plastic="DP"))                    # physics: the case's
        @test_throws ArgumentError Solver(s; grf=(; status=true))
        @test_throws ArgumentError Solver(s; bcs=(; dirichlet=fill(:fixed, 2, 2)))
        @test_throws ArgumentError Solver(s; basis=(; whch="x"))                           # an unknown option
        @test_throws ArgumentError Solver(s; solution="implicit")                          # an explicit workflow: refused
        pa, pb = Problem(s, a), Problem(s, b)                                                # one case, two solvers
        @test pa.mpts.x == pb.mpts.x && pa.basis.kind isa EP.BSplineBasis && pb.basis.kind isa EP.GimpBasis
        @test pa.time.te == 1.0 && pb.time.te == 2.0
        dp = Case(s.layout, (t = deepcopy(s.settings); t.options["material"] = Dict("plastic" => "DP"); t))
        pd = Problem(dp, a)                                                                  # one solver, two cases
        @test pd.solver === a && EP.laws(pd.mpts).plastic == "DP" && EP.laws(pa.mpts).plastic == "VM"   # each case's laws
        @test eltype(pd.mpts.s.cmp) <: EP.DruckerPrager && eltype(pa.mpts.s.cmp) <: EP.VonMises
        @test Problem(s, Solver(; dim=2, plot=(; status=false))).mpts.x == pa.mpts.x     # a solver built in code
        @test_throws ArgumentError Problem(s, Solver(; dim=3))                              # of another dimension
        @test_throws ArgumentError Solver(; material=(; plastic="DP"))                      # physics: a case's
        @test length(check_case(s).info) == 1 && length(check_case(s, a).info) == 2          # the solver's adds the time step
    end
    @testset "points: last zone wins, overrides, velocity" begin
        s = mkcase(; L=[4, 2], nel=[4, 2], ni=2,
                   zones=[Zone("all", Block([0, 0], [4, 2])), Zone("left", Block([0, 0], [1, 2])), Zone("none", Block([10, 10], [11, 11]))],
                   properties=Dict("left" => (; E=5e6, phi=0.1, v=[1.0, -2.0])))
        solver = Solver(; dim=2, plot=(; status=false))
        mesh   = EP.setup_mesh(s.layout.L, s.layout.nel, case_physics(s).bcs, solver)
        p      = case_points(s, mesh)
        left   = p.xp[1, :] .< 1
        @test size(p.xp, 2) == 32 && count(left) == 8 && all(p.zone[left] .== 2) && all(p.zone[.!left] .== 1)
        @test all(p.E[left] .== 5e6) && all(p.E[.!left] .== get_material().E) && all(p.phi[left] .== 0.1) && all(p.c0 .== get_material().c0)
        @test all(p.vp[:, left] .== [1.0, -2.0]) && all(p.vp[:, .!left] .== 0)
        x3 = case_points(mkcase(; L=[4, 1, 2], nel=[4, 1, 2], zones=s.layout.zones, properties=Dict("left" => (; v=[1.0, -2.0]))),
                         EP.setup_mesh([4.0, 1.0, 2.0], [4, 1, 2], fill(:roller, 3, 2), Solver(; dim=3, plot=(; status=false))))
        @test size(x3.xp, 2) == 64 && count(x3.zone .== 2) == 16 && all(x3.vp[:, x3.zone .== 2] .== [1.0, 0.0, -2.0])   # (x, z) velocity
    end
    @testset "TOML round trip" begin
        s = mkcase(; L=[20, 4, 10], nel=[20, 4, 10], ni=3, bcs=[:roller :roller; :roller :roller; :fixed :free], thermal=true,
                   zones=[Zone("soil", Polygon([[0, 0], [20, 0], [20, 6], [0, 6]])), Zone("rock", Block([0, 0], [20, 2])),
                          Zone("boulder", Ball([10, 4], 1.0)), Zone("lens", Polygon([[[2, 2], [8, 2], [8, 5], [2, 5]], [[4, 3], [5, 3], [5, 4]]]))],
                   properties=Dict("rock" => (; E=5e7, phi=30.0*π/180), "boulder" => (; rho0=3000.0, v=[0.0, 0.0, -1.0])),
                   material=(; c0=1.0e4, phi=25.0*π/180), solver=(; basis=(; which="gimpm", C_pf=0.37), stab=(; locking=false),
                   plot=(; what=["p"])), run=(; folder="x/y", te=1.0, workflow="elastodynamic!"))
        back = load_case(save_case(joinpath(mktempdir(), "p.toml"), s))
        @test (back.layout.L, back.layout.nel, back.layout.ni, case_physics(back).bcs, case_physics(back).thermal, Solver(back).run, case_physics(back).properties) ==
              (s.layout.L, s.layout.nel, s.layout.ni, case_physics(s).bcs, case_physics(s).thermal, Solver(s).run, case_physics(s).properties)
        z, z0 = back.layout.zones, s.layout.zones
        @test [z.name for z ∈ z] == ["soil", "rock", "boulder", "lens"] && z[3].shape.radius == 1.0 && case_physics(back).properties["boulder"].v == [0.0, 0.0, -1.0]
        @test z[4].shape.sizes == [4, 3] && z[4].shape.vertices == z0[4].shape.vertices      # several outlines
        @test case_physics(back).properties["rock"].phi == 30.0*π/180 && case_physics(back).material.phi == case_physics(s).material.phi   # radians: exact
        d = mktempdir()                                                                    # one geometry, two settings
        g = save_layout(joinpath(d, "slope.toml"), s.layout)
        dp, vm = (save_settings(joinpath(d, "$m.toml"), mkcase(; L=[20, 4, 10], nel=[20, 4, 10], ni=3, zones=s.layout.zones,
                                                           solver=(; material=(; plastic=m))).settings) for m ∈ ("DP", "VM"))
        a, b = load_case(g, dp), load_case(g, vm)
        @test (case_physics(a).laws.plastic, case_physics(b).laws.plastic) == ("DP", "VM") && points(a) == points(b)
        sb = Solver(back)
        @test case_physics(back).material.c0 == 1.0e4 && (sb.basis.which, sb.basis.C_pf, sb.stab.locking, sb.plot.what) == ("gimpm", 0.37, false, ["p"])
        @test_throws ArgumentError save_case(joinpath(mktempdir(), "r.toml"), mkcase(; L=[1, 1], nel=[1, 1], zones=[Zone("f", Region(x -> true))]))
        old = joinpath(mktempdir(), "old_case.toml")                                       # the older case format: refused
        write(old, "bcs = [[\"roller\", \"roller\"], [\"fixed\", \"free\"]]\n[domain]\nL = [1.0, 1.0]\nnel = [1, 1]\n[material]\nphi = 30.0\n")
        @test_throws ArgumentError load_case(old)
    end
    @testset "checks and time step" begin
        layered(; kw...) = mkcase(; L=[20, 10], nel=[20, 10], kw..., zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]))],
                                  properties=Dict("rock" => (; E=5e7, nu=0.2, rho0=2900.0)))
        s = layered(; solver=(; plot=(; status=false)), run=(; folder="test/spec_dt"))
        c = check_case(s)
        @test isempty(c.errors) && isempty(c.warnings) && startswith(c.info[1], "480 material points, 200 elements (2D)")
        pb = load_simulation(setup(s)).problem
        @test pb.mpts.nmp == 480
        @test estimate_dt(s)[1] ≈ EP.get_dt(pb.mpts, pb.mesh.prprt, pb.time, Inf) rtol=1e-12   # = the solver's
        ex = layered(; L=[20, 4, 10], nel=[20, 4, 10], solver=(; plot=(; status=false)), run=(; folder="test/spec_dt3"))
        pb3 = load_simulation(setup(ex)).problem
        @test startswith(check_case(ex).info[1], "$(pb3.mpts.nmp) material points, 800 elements (3D)") && pb3.mpts.nmp == 480*8
        @test estimate_dt(ex)[1] ≈ EP.get_dt(pb3.mpts, pb3.mesh.prprt, pb3.time, Inf) rtol=1e-12
        el = layered(; nel=[40, 10], solver=(; plot=(; status=false)), run=(; folder="test/spec_dt_el"))   # cells 0.5 × 1 m
        pe = load_simulation(setup(el)).problem; dt, c = estimate_dt(el)
        @test isapprox(dt, EP.get_dt(pe.mpts, pe.mesh.prprt, pe.time, Inf); rtol=1e-12) && dt ≈ 0.5*0.5/c   # the smallest cell sets it
        none = mkcase(; L=[20, 10], nel=[20, 10], zones=[Zone("away", Polygon([[30, 0], [40, 0], [40, 5]]))])
        @test length(check_case(none).errors) == 1 && length(check_case(none).warnings) == 2   # no points; outside; empty zone
        hidden = mkcase(; L=[20, 10], nel=[20, 10], zones=reverse(s.layout.zones))
        @test length(check_case(hidden).warnings) == 1                                          # rock covered by soil
        checks(c) = check_case(c, Solver(c))                                                # the numerics: the solver's checks
        @test length(checks(layered(; solver=(; material=(; elastic="hypoelastic")), run=(; workflow="elastoquasistatic!"))).errors) == 1
        @test length(checks(layered(; solver=(; solution="explicit"), run=(; workflow="elastoquasistatic!"))).errors) == 1
        @test_throws ArgumentError Solver(layered(; solver=(; solution="implicit")))       # an explicit workflow: refused
        @test length(checks(layered(; solver=(; dtype=(; precision="32-bit precision (or single precision)"), basis=(; which="gimpm")))).warnings) == 1
        @test length(checks(layered(; run=(; tg=20.0))).warnings) == 1
        @test length(checks(mkcase(; L=[20, 10], nel=[10, 10], zones=s.layout.zones, run=(; workflow="elastoquasistatic!"))).warnings) == 1
        @test isempty(check_case(layered(; run=(; tg=20.0))).warnings)                      # a case alone: its physics
    end
    @testset "zone constants, boundaries and options reach the problem" begin
        s = mkcase(; L=[20, 10], nel=[20, 10], bcs=[:roller :roller; :fixed :free],
                   zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]))],
                   properties=Dict("rock" => (; E=5e7, nu=0.2, rho0=2900.0)),
                   solver=(; basis=(; which="gimpm", C_pf=0.37, domain="undeformed"), stab=(; locking=false), material=(; plastic="VM"), plot=(; status=false)),
                   run=(; folder="test/case_zones"))
        pb = load_simulation(setup(s)).problem; solver = pb.solver
        m = pb.mpts
        @test solver.basis == (; which="gimpm", trsfr=get_default().basis.trsfr, C_pf=0.37, domain="undeformed") && !solver.stab.locking &&
              solver.stab.damping == get_default().stab.damping && EP.laws(m).plastic == "VM"
        rock = getindex.(m.x, 2) .< 2.0
        K, G = EP.get_elastic_moduli(5e7, 0.2); K0, G0 = EP.get_elastic_moduli(get_material().E, get_material().nu)
        @test all(c -> (c.Kc, c.Gc) == (K, G), m.s.cmp[rock]) && all(m.s.ρ₀[rock] .== 2900.0)
        @test all(c -> (c.Kc, c.Gc) == (K0, G0), m.s.cmp[.!rock]) && all(m.s.ρ₀[.!rock] .== get_material().rho0)
        base = [i for i ∈ eachindex(pb.mesh.x) if pb.mesh.x[i][2] == 0]
        @test all(pb.mesh.s.bcs.status[:, base])
        F, ℓ₀, ℓ = EP.SMatrix{2,2}(2.0, 0.0, 0.0, 0.5), EP.SVector(0.25, 0.25), EP.SVector(0.3, 0.2)   # GIMP's domain updates
        @test pb.basis.kind.domain isa EP.UndeformedDomain && EP.particle_domain(pb.basis.kind, F, ℓ₀, ℓ) == ℓ₀
        @test EP.particle_domain(EP.GimpBasis{Int,2}(EP.StretchDomain()), F, ℓ₀, ℓ) ≈ [0.5, 0.125]          # ℓ₀ ∘ diag(U)
        @test EP.particle_domain(EP.LinearBasis{Int,2}(), F, ℓ₀, ℓ) === ℓ                                  # other kinds: kept
        @test_throws ArgumentError EP.get_domain("cpdi")
    end
    @testset "free option values" begin
        @test EP.option_value(["basis", "C_pf"], 0.37) == 0.37 && EP.option_value(["grf", "param", "Nₕ"], "2000") === 2000
        @test EP.option_value(["grf", "param", "Iₓ"], "3, 4, 5") == [3.0, 4.0, 5.0] && EP.option_value(["stab", "locking"], "false") === false
        @test_throws ArgumentError EP.option_value(["grf", "param", "Nₕ"], 2.5)                    # an integer option
        @test_throws ArgumentError EP.option_value(["basis", "which"], "foo")                      # not a numeric option
        @test EP.option_value(["plot", "what"], ["p", "plastic strain"]) == ["p", "epII"]          # a key, or an older case's name
    end
    @testset "simulation: build, save, load, execute!" begin
        dir = joinpath(mktempdir(), "mycase")
        c = mkcase(; L=[4, 2], nel=[4, 2], zones=[Zone("a", Block([0, 0], [4, 1]))], solver=(; plot=(; status=false)),
                   run=(; folder=dir * "/", te=0.2, tg=0.1, tep=0.0, workflow="elastodynamic!"))
        sim = Simulation(Problem(c, Solver(c)); folder=EP.case_folder(c))
        @test sim.folder == dir && !isdir(dir) && fieldnames(Simulation) == (:problem, :folder)   # built, nothing written
        @test Simulation(sim.problem; folder="a/b").folder == joinpath(EP.self.sys.dump, "a", "b")   # a name: under dump/
        @test_throws ArgumentError Solver(c; run=(; workflow="nope"))
        @test_throws ArgumentError Solver(c; run=(; lstps=0))
        jld2 = setup(c)
        @test jld2 == joinpath(dir, "mycase_simulation.jld2") && isfile(joinpath(dir, "mycase_case.toml"))
        l = load_simulation(jld2); @test (l.problem.solver.run, l.folder) == (sim.problem.solver.run, sim.folder)
        saved = jldopen(f -> f["simulation"].problem.solver, jld2)
        @test saved.cairn === nothing && keys(saved.backend) == (:select,)                  # no kernels or devices
        @test l.problem.solver.cairn isa NamedTuple && l.problem.solver.backend.exec.backend isa EP.CPU   # rebuilt
        state(p) = (p.mpts.x, getfield.(p.mpts.s.σᵢ, :p), p.mpts.s.v)
        s0 = state(sim.problem); r = execute(jld2)
        @test state(r.problem) != s0 && state(load_simulation(jld2).problem) == s0               # ran a copy, file untouched
        @test r.problem.time.t[1] == 0.2                                                          # to te, though te < plot.freq
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
        bc = EP.block([8, 4]; plot=(; status=false, freq=1.0, dpi=100, what=[]))     # quasi-static: implicit
        @test Solver(bc).solution == "implicit" && Solver(bc).run.workflow == "elastoquasistatic!"
        @test_throws ArgumentError elastodynamic!(Problem(bc, Solver(bc)))                   # explicit: refuses it
        @test_throws ArgumentError Solver(bc; solution="implicit", run=(; workflow="elastodynamic!"))   # an explicit one: refused
        qs = Simulation(Problem(bc, Solver(bc; run=(; lstps=3))))
        @test length(only(case_workflows(qs))(qs.problem)) == 3                                   # 3 load steps
        @test_throws ArgumentError Solver(c; run=(; folder="elsewhere"))                    # the folder: the case's
        gimp  = Solver(c; basis=(; which="gimpm"), run=(; te=0.1))
        other = Case(c.layout, (t = deepcopy(c.settings); t.run["folder"] = dir * "_other"; t))
        j1, j2 = setup(c, gimp), setup(other, gimp)                    # one solver, two cases
        @test (dirname(j1), dirname(j2)) == (dir, dir * "_other")                                 # each in its own folder
        @test EP.solver_config(Solver(load_case(joinpath(dir, "mycase_case.toml")))) == EP.solver_config(gimp)   # as solved
    end
    @testset "containers adapt: host arrays and back, unchanged" begin
        c = EP.slump([64.1584, 16.0396], [16, 4]; plot=(; status=false), basis=(; which="gimpm", trsfr="apic"))
        p = Problem(c, Solver(c))
        for x ∈ (p.mpts, p.mesh, p.basis)                                                  # the array types are parameters
            y = EP.Adapt.adapt(Array, x)
            @test typeof(y) == typeof(x)
            @test all(f -> getfield(x, f) isa EP.AbstractTransfer || getfield(y, f) == getfield(x, f), fieldnames(typeof(x)))
        end
        @test EP.Adapt.adapt(Array, p.basis).transfer.Bᵢⱼ == p.basis.transfer.Bᵢⱼ
        @test isbits(p.mesh.prprt) && p.basis.e2e isa EP.Adjacency                         # plain value; index vectors
    end
    @testset "workflows and a short run" begin
        s = mkcase(; L=[20, 4, 10], nel=[20, 4, 10], zones=[Zone("soil", Block([0, 0], [20, 6]))],
                   solver=(; plot=(; status=false)), run=(; folder="test/case_run", te=0.2, tg=0.1, tep=0.0, workflow="elastodynamic!"))
        @test case_workflows(Solver(s).run.workflow, Solver(s).run.lstps) == Function[elastodynamic!]
        @test case_workflows(EP.get_run().workflow, 1) == Function[elastodynamic!, elastoplastic!]
        @test length(case_workflows("elastoquasistatic! (plastic)", 10)) == 1
        @test isfile(execute!(setup(s)))                                               # 3D, its own times
    end
end
