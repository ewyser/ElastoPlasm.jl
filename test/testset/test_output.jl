@testset "+ $(basename(@__FILE__))" verbose = true begin
    # the HDF5 time series of a run (output.jl, output_plots.jl), on a small layered case
    import HDF5
    EP   = ElastoPlasm
    out(; kw...) = (; status=true, every=0.5, kw...)
    case(name; output=out(), te=1.0, tep=0.5, workflow="elastodynamic! + elastoplastic!", solver=(;)) =
        Case(; L=[20, 10], nel=[20, 10], bcs=[:roller :roller; :fixed :free],
               zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]); E=5e7)],
               solver=(; plot=(; status=false), output, solver...), run=(; folder="test/$name", te, tg=te, tep, workflow))
    output_of(j) = replace(j, "_simulation.jld2" => "_output.h5")

    @testset "configuration" begin
        s = get_solver(; output=(; status=true, what=(; point=(; solid=["p"]))))
        @test s.output.what.point.solid == ["p"] && s.output.what.point.base == get_default().output.what.point.base   # per-field merge
        @test_throws ArgumentError EP.option_value(["output", "what", "point", "solid"], ["pp"])
        c = case("out_cfg"; output=out(what=(; point=(; solid=["p", "q"]))))
        @test load_case(save_case(joinpath(mktempdir(), "c.toml"), c)).solver.output.what.point.solid == ["p", "q"]
    end
    @testset "a run recorded" begin
        j = setup_simulation(case("out_run"))
        execute!(j; workflows=[elastodynamic!])
        f = output_of(j)
        @test output_series(f, "point/solid/p"; point=1).t == [0.0, 0.5, 1.0]
        execute!(j; workflows=[elastoplastic!])                                    # continues the saved state: appended
        p = load_simulation(j).problem
        HDF5.h5open(f, "r") do h                                                   # plain HDF5, no ElastoPlasm types
            t = read(h["frames/t"])
            @test t == [0.0, 0.5, 1.0, 1.5] && issorted(t) && read(h["frames/loadstep"]) == zeros(UInt8, 4)
            @test size(h["point/solid/u"]) == (2, p.mpts.nmp, 4) && size(h["point/solid/p"]) == (p.mpts.nmp, 4) &&
                  size(h["point/solid/dev"]) == (6, p.mpts.nmp, 4) && eltype(h["point/solid/p"]) == Float64
            @test h["point/base/x"][:, :, end] == reduce(hcat, p.mpts.x) && h["point/solid/p"][:, end] == getfield.(p.mpts.s.σᵢ, :p)   # last frame = saved state
            @test read(HDF5.attributes(h)["nmp"]) == p.mpts.nmp && occursin("[[zone]]", read(HDF5.attributes(h)["case"]))
            a = HDF5.attributes(h["point/solid/p"])                                # each field says what it is
            @test (read(a["name"]), read(a["unit"]), read(a["label"])) == ("pressure", "Pa", "\$p\$") &&
                  read(HDF5.attributes(h["frames/t"])["unit"]) == "s" && !haskey(HDF5.attributes(h["point/solid/u"]), "label")
        end
        s = output_series(f, "point/solid/p"; point=7)
        @test s.t == [0.0, 0.5, 1.0, 1.5] && s.values[end] == p.mpts.s.σᵢ[7].p && !s.loadstep
        @test output_series(f, "point/solid/u"; point=7, component=2).values[end] == p.mpts.s.u[7][2]
        fr = output_frame(f, 3)
        @test fr.t == 1.0 && fr.workflow == "elastodynamic!" && fr.fields["point/solid/epII"] == HDF5.h5open(h -> h["point/solid/epII"][:, 3], f)
        @test load_output(f).frames.workflow == ["initial", "elastodynamic!", "elastodynamic!", "elastoplastic!"]
        @test load_output(f).info["point/solid/u"] == (; name="displacement", unit="m", label=nothing)
        execute(j; workflows=[elastodynamic!])                                     # a copy of the saved state: a new file
        @test output_series(f, "point/solid/p"; point=1).t[1] == load_simulation(j).problem.time.t[1]
    end
    @testset "selection, phases" begin
        c = case("out_sel"; output=out(what=(; point=(; solid=["p"], thermal=["T"]), mesh=(; solid=["v"]))))
        j = setup_simulation(c)
        @test_logs (:warn, r"skipped fields") match_mode=:any execute!(j; workflows=[elastodynamic!])
        HDF5.h5open(output_of(j), "r") do h
            @test !haskey(h, "point/thermal") && keys(h["point/solid"]) == ["p"]
            @test size(h["mesh/solid/v"]) == (2, load_simulation(j).problem.mesh.prprt.nno[end], length(h["frames/t"]))
        end
        p = execute!(Simulation(Problem(case("out_mem"))); workflows=[elastodynamic!], output=(f = joinpath(mktempdir(), "o.h5"))).problem   # in memory
        @test isfile(f) && output_series(f, "point/solid/p"; point=1).t[end] == p.time.t[1]
    end
    @testset "load steps" begin
        b = setup_simulation(EP.block_case([16, 9]; plot=(; status=false), output=out(), folder="test/out_block"))
        execute!(b; workflows=[(a...) -> elastoquasistatic!(a...; lstps=2, G=10.0)])
        s = output_series(output_of(b), "point/solid/p"; point=1)
        @test s.t == [0.0, 1.0, 2.0] && s.loadstep
    end
    @testset "visualise after the run" begin
        j = setup_simulation(case("out_plot"))
        execute!(j; workflows=[elastodynamic!, elastoplastic!])
        f, d = output_of(j), mktempdir()
        @test plot_output(f, "point/solid/p"; frame=2, file_out=joinpath(d, "p.png")) isa EP.Plots.Plot && isfile(joinpath(d, "p.png"))
        @test plot_output(f, "point/solid/u") isa EP.Plots.Plot                     # a vector: its magnitude
        @test_throws ArgumentError plot_output(f, "point/solid/dev")                # 6 components: say which
        @test isfile(animate_output(f, "point/solid/p"; fps=2, file_out=joinpath(d, "p.gif")))
        @test plot_series(f, "point/solid/p"; point=5) isa EP.Plots.Plot
    end
end
