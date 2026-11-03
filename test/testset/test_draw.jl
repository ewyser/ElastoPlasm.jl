@testset "+ $(basename(@__FILE__))" verbose = true begin
    # a case written in code: the layout, then its settings in ElastoPlasm's terms, read at once (its errors here)
    function mkcase(; L, nel, ni=2, zones=Zone[], kwargs...)
        c = Case(Layout(; L, nel, ni, zones), case_settings(; kwargs...))
        case_physics(c)
        return c
    end
    # ElastoPlasm's side of the drawing window, which is EasyDraw's (lib/EasyDraw, tested by test_easydraw.jl): its
    # catalog, a case's layout in the window and back, and a case drawn through the editor, driven as a user would (its options,
    # run settings, checks and files); then draw_case's window, offscreen
    EP, ED = ElastoPlasm, ElastoPlasm.EasyDraw
    editor(case) = ED.Editor(case.layout, case.settings, get_catalog(); EP.case_hooks()...)
    click!(e, x, z) = (ED.move!(e, x, z); ED.press!(e, x, z); ED.release!(e))
    opts(c) = EP.config_dict(EP.case_options(c.settings.options))                          # a case's options, as values

    @testset "the catalog: every solver option, the zone constants, the run" begin
        d = get_catalog()                                                                   # plain data, TOML-shaped:
        @test d isa Dict && occursin("[[options]]", sprint(EP.TOML.print, d))
        c = ED.Catalog(d)                                                                   # as EasyDraw reads it
        leaves(o, path=String[]) = o isa NamedTuple ? reduce(vcat, (leaves(getfield(o, k), [path; string(k)]) for k ∈ keys(o)); init=String[]) :
                                   isempty(o[2]) ? String[] : [join(path, ".")]
        @test sort([f.key for f ∈ c.options]) == sort(setdiff(leaves(get_option()), ["dtype.T0", "dtype.bits", "solution"]))
        @test all(f -> f.default == EP.default_at(split(f.key, ".")), c.options)          # get_default()'s values
        @test first(ED.sections(c)) == "basis" && ED.field(c.options, "grf.param.Nₕ").kind == :integer
        @test ED.field(c.options, "grf.param.Iₓ").kind == :numeric && ED.field(c.options, "basis.which").kind == :menu
        w = ED.field(c.options, "plot.what"); @test w.kind == :list && w.names[findfirst(==("epII"), w.values)] == "plastic strain"
        @test [f.key for f ∈ c.zone] == ["c0", "cr", "phi", "T", "E", "nu", "rho0"] && all(f -> f.default == getfield(get_material(), Symbol(f.key)), c.zone)
        @test ED.field(c.zone, "phi").show_if == Dict("material.plastic" => ["DP"]) && c.side == "roller"
        @test [f.key for f ∈ c.run] == ["workflow", "te", "tg", "tep", "lstps"] && ED.field(c.run, "workflow").values == collect(EP.get_workflows())
        back = ED.read_catalog(ED.write_catalog(joinpath(mktempdir(), "c.toml"), d))                # as TOML
        @test all(k -> ED.field_dict.(getfield(back, k)) == ED.field_dict.(getfield(c, k)), (:zone, :options, :run))
    end
    @testset "a case's layout in the window and back, unchanged" begin
        case = mkcase(; L=[20, 4, 10], nel=[20, 4, 10], ni=3, bcs=[:roller :roller; :free :fixed; :fixed :free], thermal=true,
                      zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]))],
                      properties=Dict("rock" => (; E=5e7, v=[0.0, -1.0])),
                      material=(; c0=1.0e4), solver=(; basis=(; which="gimpm"), plot=(; what=["p", "epII"])), run=(; folder="x/y", te=1.0))
        e = editor(case); back = Case(ED.current_layout(e), ED.current_settings(e))         # opened, nothing edited
        @test (back.layout.L, back.layout.nel, back.layout.ni, case_physics(back).bcs, case_physics(back).thermal, case_physics(back).material, Solver(back).run, case_physics(back).properties) ==
              (case.layout.L, case.layout.nel, case.layout.ni, case_physics(case).bcs, case_physics(case).thermal, case_physics(case).material, Solver(case).run, case_physics(case).properties)   # y sides carried
        @test opts(back) == opts(case)                                                      # (a Dict's order aside)
        @test [z.name for z ∈ back.layout.zones] == [z.name for z ∈ case.layout.zones]
        @test points(back) == points(case)                                                  # a block comes back as its rectangle
    end
    @testset "the window's points are the simulation's, to the bit" begin
        for c ∈ (slump_case([64.1584, 16.0396], [40, 10]), EP.collapse_case([64, 36]), EP.collision_case([1.0, 1.0], [20, 20]),
                 mkcase(; L=[10.3, 7.1], nel=[13, 9], ni=3, zones=[Zone("ball", Ball([5.15, 3.55], 2.0))]),
                 mkcase(; L=[8, 4], nel=[8, 4], zones=[Zone("45°", Polygon([[0, 0], [4, 0], [0, 4]]))]))   # an edge through centres
            solver = Solver(; dim=2, plot=(; status=false))
            p = case_points(c, EP.setup_mesh(c.layout.L, c.layout.nel, case_physics(c).bcs, solver))
            sim = sort(collect(zip(p.xp[1, :], p.xp[2, :], p.zone)))
            q = points(c); x, z, r = ED.preview_points(editor(c))                          # EasyDraw's; the window's preview
            @test sim == sort(collect(zip(q.x, q.z, q.zone))) && sim == sort(collect(zip(x, z, r)))
            @test count(>(0), EP.case_zones(c)) == size(p.xp, 2)                           # check_case counts the same
        end
    end
    @testset "layout files as cases, cases in EasyDraw's functions" begin
        case = mkcase(; L=[20, 4, 10], nel=[20, 4, 10], ni=3, bcs=[:roller :roller; :roller :roller; :fixed :free], thermal=true,
                      zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]))],
                      properties=Dict("soil" => (; phi=0.5), "rock" => (; E=5e7, v=[0.0, -1.0])),
                      material=(; c0=1.0e4), solver=(; basis=(; which="gimpm")), run=(; te=1.0))
        d = mktempdir(); f = save_case(joinpath(d, "x_layout.toml"), case)              # for EasyDraw alone
        s = ED.load_settings(f)
        @test s.zones["soil"]["phi"] == 0.5 && s.defaults["c0"] == 1.0e4 && s.sides["bottom"] == "fixed"   # SI
        back = load_case(f)                                                                # and back, as a case
        @test (back.layout.L, back.layout.nel, back.layout.ni, case_physics(back).bcs, case_physics(back).thermal, case_physics(back).material, Solver(back).run, case_physics(back).properties) ==
              (case.layout.L, case.layout.nel, case.layout.ni, case_physics(case).bcs, case_physics(case).thermal, case_physics(case).material, Solver(case).run, case_physics(case).properties)
        @test opts(back) == opts(case)
        @test_logs ED.load_layout(save_case(joinpath(d, "x_case.toml"), case))              # a case file is a layout file
        e = ED.Editor(ED.Layout(; L=[10, 5], nel=[10, 5]), ED.Settings(), get_catalog())   # drawn with EasyDraw alone, ElastoPlasm's catalog
        ED.new_zone!(e, :rectangle); click!(e, 0.0, 0.0); click!(e, 10.0, 2.0)
        ED.set_prop!(e, "phi", "30"); ED.set_bc!(e, "bottom", "fixed"); ED.set_option!(e, "basis.which", 2)
        @test ED.create!(e, joinpath(d, "drawn"))                                           # EasyDraw's Create…: a layout file
        c = load_case(joinpath(d, "drawn", "drawn_layout.toml"))                            # a case
        @test case_physics(c).properties["zone 1"].phi ≈ 30π/180 && case_physics(c).bcs[2, 1] == :fixed && Solver(c).basis.which == "gimpm"
        @test EP.case_folder(c) == joinpath(d, "drawn")                                     # no run folder: the file's
        @test isempty(check_case(c).errors)
        @test Problem(c, Solver(c)).mpts.nmp == 80                                                     # 10 × 2 elements, 2 × 2 points each
        @test zone_at(case, (10.0, 1.0)) == 2 && zone_at(case, (10.0, 3.0, 1.0)) == 2 && zone_at(case, (10.0, 8.0)) == 0   # (x, y, z) too
        xs, zs = 0.5:1.0:19.5, 0.5:1.0:9.5
        @test zone_index(case, xs, zs)[:, 1] == fill(2, 20) && zone_index(case, xs, zs)[:, 10] == fill(0, 20)
        E, c0 = property(case, "E", xs, zs; default=get_material().E), property(case, "c0", xs, zs)
        @test E[1, 1] == 5e7 && E[1, 4] == get_material().E && all(==(1.0e4), c0)                # the case's material where not set
        @test_throws ArgumentError property(case, "v", xs, zs)                             # not a number
        @test length(points(case).x)*4*3 == count(>(0), EP.case_zones(case))               # the x–z plane; the case counts 3D
        @test sort(unique(r.zone for r ∈ partition(case).regions)) == [0, 1, 2]
        @test occursin("Physical Surface(\"rock\")", read(write_geo(joinpath(d, "x.geo"), case), String))
        @test_throws ArgumentError partition(mkcase(; L=[1, 1, 1], nel=[1, 1, 1], zones=[Zone("b", Ball([0.5, 0.5, 0.5], 0.2))]))   # a 3D shape
    end
    @testset "option values: a case file and the window read them alike" begin
        for f ∈ ED.Catalog(get_catalog()).options                                          # ElastoPlasm's option_value, EasyDraw's field_value
            xs = f.kind == :list ? Any[f.default, f.values[1:1], string.(f.values), f.names, ["zzz"], "p"] :
                 Any[f.default, f.values..., string.(f.values)..., "abc", "2.5", "3", "1, 2, 3", "Inf", "NaN", 2.5, 3, true]
            for x ∈ xs
                a = try EP.option_value(split(f.key, "."), x) catch; :error end
                b = try ED.field_value(f, x) catch; :error end
                isequal(a, b) || @info "differ" f.key x a b
                @test isequal(a, b)
            end
        end
    end
    @testset "a case drawn: zone constants, options, run, checks, files" begin
        e = editor(mkcase(; L=[13.0, 7.0], nel=[13, 7], bcs=[:roller :roller; :fixed :free]))
        @test e.sides == Dict("left" => "roller", "right" => "roller", "bottom" => "fixed", "top" => "free")
        ED.set_bc!(e, "top", "fixed"); ED.set_domain_field!(e, "width", "4")
        ED.new_zone!(e); foreach(p -> click!(e, p...), [(0, 0), (13, 0), (13, 3), (0, 3)]); ED.key!(e, "enter")
        ED.set_prop!(e, "E", "50"); ED.set_prop!(e, "nu", "0.2"); ED.set_prop!(e, "phi", "30")
        @test e.regions[1].props == Dict("E" => 5e7, "nu" => 0.2, "phi" => 30π/180)
        @test occursin("default 20.0", ED.panel_json(e))                                    # c₀'s default, get_material()'s, in kPa
        ED.set_option!(e, "basis.which", findfirst(==("gimpm"), ED.field(e.catalog.options, "basis.which").values))
        @test e.options["basis"]["which"] == "gimpm"
        ED.type_option!(e, "grf.param.σ", "7500"); @test e.options["grf"]["param"]["σ"] == 7500.0
        ED.tick_value!(e, "plot.what", findfirst(==("epII"), ED.field(e.catalog.options, "plot.what").values), true)
        @test e.options["plot"]["what"] == ["p", "epII"]
        e.section = "plot"; @test occursin("\"plastic strain\"", ED.panel_json(e))           # shown by its readable name
        ED.set_option!(e, "material.plastic", 2); @test !occursin("φ [°]", ED.panel_json(e))  # von Mises: no friction angle
        ED.pick_run!(e, "workflow", findfirst(==("elastoquasistatic!"), collect(EP.get_workflows())))
        ED.set_run!(e, "lstps", "20"); ED.set_run!(e, "lstps", "-1")
        @test e.run["workflow"] == "elastoquasistatic!" && e.run["lstps"] == 20
        @test occursin("20 load steps", join(ED.checks(e).info))                              # check_case's
        ED.setoption!(e.options, "dtype", ["precision"], "32-bit precision (or single precision)")   # a warning: Create… twice
        @test !ED.ready!(e) && occursin("again", e.status) && ED.ready!(e)
        folder = joinpath(mktempdir(), "layers")
        @test ED.create!(e, folder * "/") && e.finished && isfile(joinpath(folder, "layers_case.toml"))   # made, the case in it
        c = Case(e.result.layout, e.result.settings)
        sc = Solver(c)
        @test EP.case_folder(c) == folder && c.layout.L == [13.0, 4.0, 7.0] && sc.run.lstps == 20 && sc.basis.which == "gimpm"   # its runs go there
        @test case_physics(c).bcs == [:roller :roller; :roller :roller; :fixed :fixed] && case_physics(c).properties["zone 1"].E == 5e7
        @test Solver(load_case(joinpath(folder, "layers_case.toml"))).run == sc.run          # the file is the case
        e2 = editor(mkcase(; L=[1, 1], nel=[2, 2]))
        ED.load!(e2, joinpath(folder, "layers_case.toml"))                                  # Start from a (case) file
        @test length(e2.regions) == 1 && e2.width == 4.0 && e2.run["lstps"] == 20 && e2.folder == folder
        @test opts(Case(ED.current_layout(e2), ED.current_settings(e2))) == opts(c)
    end
    ENV["QT_QPA_PLATFORM"] = get(ENV, "QT_QPA_PLATFORM", "offscreen"); ENV["QT_QUICK_BACKEND"] = "software"
    @testset "draw_case's window, offscreen" begin
        d = mktempdir()
        case = mkcase(; L=[20.0, 10.0], nel=[20, 10], bcs=[:roller :roller; :fixed :free],
                     zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]))], properties=Dict("rock" => (; E=5e7)))
        @test isnothing(draw_case(case; shots=(d, [1, 2, 3, 4])))                               # closes itself; no Create…: nothing
        @test all(k -> isfile(joinpath(d, "step$k.png")), 1:4)
    end
end
