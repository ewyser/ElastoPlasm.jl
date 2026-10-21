@testset "+ $(basename(@__FILE__))" verbose = true begin
    # the drawing window: its Qt-free part (ext/DrawQMLExt/editing.jl: the Editor and the JSON the view reads), loaded
    # on its own and driven as a user would, everywhere; then the QML window itself, offscreen, where QML is installed.
    EP = ElastoPlasm
    ed = @eval module DrawEditing
        using ElastoPlasm
        include(joinpath(pkgdir(ElastoPlasm), "ext", "DrawQMLExt", "editing.jl"))
    end
    L1, L2 = 64.1584, 64.1584/4
    click!(e, x, z) = (ed.move!(e, x, z); ed.press!(e, x, z); ed.release!(e))

    @testset "drawing a zone" begin
        e = ed.Editor(Case(; L=[L1, L2], nel=[40, 10]))
        ed.new_zone!(e)
        foreach(p -> click!(e, p...), [(0.1, 0.05), (64.1, 0.1), (64.0, 1.95), (22.5, 12.05), (0.05, 12.0), (0.2, 0.1)])   # last: closes
        ed.set_prop!(e, "phi", "7.5")
        @test ed.done!(e) && e.result isa Case && length(e.result.zones) == 1
        v = e.result.zones[1].shape.vertices
        @test length(v) == 5 && v[1] == [0.0, 0.0] && v[4][2] ≈ 15*L2/20                  # snapped to the sub-cell lattice
        @test e.result.zones[1].props[:phi] ≈ 7.5*π/180 rtol=1e-14
    end
    @testset "closing on the first vertex" begin
        e = ed.Editor(Case(; L=[L1, L2], nel=[40, 10]))
        ed.new_zone!(e); click!(e, 10.0, 2.0); click!(e, 30.0, 2.0)
        first = copy(e.drawing[1]); near = first .+ 0.3                                     # within the snap radius
        ed.move!(e, near...); @test e.cursor[2] != :close                                   # 2 vertices: no closing yet
        ed.press!(e, near...); ed.release!(e); @test length(e.drawing) == 2                 # and no repeated first vertex
        click!(e, 30.0, 8.0)
        ed.move!(e, near...)
        @test e.cursor == (first, :close) && occursin("close", e.readout)                  # snapped onto it
        @test occursin("\"r\":9", ed.overlay_json(e, 800.0, 600.0))                           # drawn as the closing marker
        ed.press!(e, near...); ed.release!(e)
        @test isnothing(e.drawing) && length(e.regions) == 1                                # closed, exactly on it:
        @test e.regions[1].vertices[1] == first && length(e.regions[1].vertices) == 3       # no extra vertex
    end
    @testset "snap radius follows the zoom; the cursor's colour and shape" begin
        e = ed.Editor(Case(; L=[L1, L2], nel=[40, 10], zones=[Zone("base", Block([0, 0], [L1, 4]))]))
        v = [0.0, 4.0]; @test v ∈ e.regions[1].vertices
        r0 = ed.tol(e); ed.zoom!(e, 5.0, 4.0, 0.1); @test ed.tol(e) ≈ 0.1r0               # 10× in: a 10× smaller radius
        ed.move!(e, (v .+ 0.5r0)...);  @test e.cursor[2] == :grid                          # no longer on the vertex
        ed.move!(e, (v .+ 0.05r0)...); @test e.cursor[2] == :vertex
        ed.fit!(e)
        cur() = ed.cursor_json(e, (x, z) -> [x, z])
        ed.move!(e, v...);      @test cur().color == "#9ca3af" && cur().shape == "square"   # not drawing: grey
        ed.move!(e, 30.0, 4.0); @test cur().color == "#9ca3af" && cur().shape == "diamond"  # on the base's top edge
        ed.new_zone!(e); ed.move!(e, 30.0, 9.0)
        @test cur().color == "#dc2626" && cur().shape == "ring"                             # drawing: red
        @test occursin("\"drawing\":true", ed.panel_json(e))                              # the crosshair, the highlighted button
    end
    @testset "the picked vertex, typed" begin
        e = ed.Editor(Case(; L=[20.0, 10.0], nel=[20, 10],
                           zones=[Zone("A", Block([0, 0], [20, 4])), Zone("B", Polygon([[0, 4], [20, 4], [20, 8], [0, 8]]))]))
        A, B = () -> e.regions[1].vertices, () -> e.regions[2].vertices
        @test isnothing(e.picked) && occursin("\"vertex\":null", ed.panel_json(e))         # nothing picked: no row
        ed.move!(e, 20.0, 4.0); ed.press!(e, 20.0, 4.0); ed.release!(e)                     # click the shared corner
        @test e.picked == [20.0, 4.0] && e.selected == 2 && occursin("\"vertex\":{\"x\":\"20.0\",\"z\":\"4.0\"}", ed.panel_json(e))
        ed.set_picked!(e, "z", "5")
        @test [20.0, 5.0] ∈ A() && [20.0, 5.0] ∈ B() && [20.0, 4.0] ∉ A() && e.picked == [20.0, 5.0]   # moved in both zones
        ed.set_picked!(e, "x", "abc"); @test e.picked == [20.0, 5.0]                         # not a number: ignored
        ed.undo!(e); @test [20.0, 4.0] ∈ B() && occursin("\"vertex\":null", ed.panel_json(e))   # the picked one is gone
        ed.move!(e, 0.0, 8.0); ed.press!(e, 0.0, 8.0); ed.release!(e); n = length(B())
        ed.remove_picked!(e); @test length(B()) == n - 1 && isnothing(e.picked)              # deleted
        @test occursin("\"picked\":null", ed.overlay_json(e, 800.0, 600.0))
    end
    @testset "new zones: rectangle, triangle" begin
        e = ed.Editor(Case(; L=[20.0, 10.0], nel=[20, 10], zones=[Zone("base", Block([0, 0], [20, 2]))]))
        ed.new_zone!(e, :rectangle); click!(e, 4.0, 2.0); click!(e, 4.0, 6.0)               # not opposite corners
        @test length(e.regions) == 1 && occursin("opposite", e.status)
        click!(e, 10.0, 6.0)
        @test length(e.regions) == 2 && isnothing(e.drawing) && Set(e.regions[2].vertices) == Set([[4.0, 2.0], [10.0, 2.0], [10.0, 6.0], [4.0, 6.0]])
        @test [4.0, 2.0] ∈ e.regions[1].vertices && [10.0, 2.0] ∈ e.regions[1].vertices     # corners on the base: shared
        ed.new_zone!(e, :triangle); click!(e, 12.0, 2.0); click!(e, 18.0, 2.0)
        ed.move!(e, 15.0, 7.0); @test length(ed.sketch_px(e, (x, z) -> [x, z])) == 4          # closed back to its first corner
        click!(e, 15.0, 7.0)
        @test length(e.regions) == 3 && isnothing(e.drawing) && length(e.regions[3].vertices) == 3   # closes on its third
        @test_throws ArgumentError ed.new_zone!(e, :hexagon)
    end
    @testset "circles (Ball zones)" begin
        case = EP.collision_case([1.0, 1.0], [20, 20])                                      # two disks, with velocities
        e = ed.Editor(case)
        @test all(ed.iscircle, e.regions) && e.regions[1].centre == [0.25, 0.25] && e.regions[1].radius == 0.1
        x, z, r = ed.preview_points(e); mp = Problem(case).mpts
        a, b = sort([[x[i], z[i]] for i ∈ eachindex(x)]), sort([[p[1], p[2]] for p ∈ mp.x])
        @test length(a) == mp.nmp && all(isapprox.(a, b; atol=1e-12))                       # the same points
        back = ed.current_case(e).zones
        @test all(i -> back[i].shape isa Ball && back[i].shape.centre == case.zones[i].shape.centre &&
                       back[i].shape.radius == case.zones[i].shape.radius && back[i].props[:v] == case.zones[i].props[:v], 1:2)
        ed.new_zone!(e, :circle); click!(e, 0.5, 0.5); click!(e, 0.6, 0.5)                          # centre, then a point on it
        c = e.regions[3]; @test ed.iscircle(c) && c.centre == [0.5, 0.5] && c.radius ≈ 0.1 && e.selected == 3
        ed.move!(e, 0.5, 0.5); ed.press!(e, 0.5, 0.5); ed.move!(e, 0.7, 0.6); ed.release!(e)  # drag the centre
        @test c.centre ≈ [0.7, 0.6] && c.radius ≈ 0.1
        h = ed.radius_handle(c); ed.move!(e, h...); ed.press!(e, h...); ed.move!(e, 0.95, 0.6); ed.release!(e)
        @test c.radius ≈ 0.25                                                               # the handle sets the radius
        ed.set_circle!(e, "r", "0.15"); ed.set_circle!(e, "r", "-1"); @test c.radius == 0.15   # typed; ≤ 0 ignored
        ed.undo!(e); @test e.regions[3].radius ≈ 0.25
        @test occursin("\"circle\":{\"x\":", ed.panel_json(e)) && occursin("\"handles\"", ed.overlay_json(e, 800.0, 600.0))
    end
    @testset "a picture to draw on: its scale, and the domain fitted to it" begin
        e = ed.Editor(Case(; L=[30.0, 15.0], nel=[30, 15]); step=1)                      # h = 1 m
        ed.set_backdrop!(e, "section.png")
        @test occursin("\"backdrop\":null", ed.background_json(e, 800.0, 600.0))           # drawn once its size is known
        ed.backdrop_size!(e, 800, 400); b = e.backdrop
        @test (b.x0, b.z0, b.width, ed.height(b)) == (0.0, 0.0, 30.0, 15.0)                 # until scaled: the domain's width
        @test occursin("\"backdrop\":{\"url\":\"file://", ed.background_json(e, 800.0, 600.0))
        ed.set_backdrop_field!(e, "width", "24"); @test ed.height(b) == 12.0 && e.L == [24.0, 12.0] && e.nel == [24, 12]
        # a scale bar drawn 5 m long on the picture, as it shows at width 24: its ends 6 m apart; true width 20 m
        ed.scale_bar!(e); ed.key!(e, "escape"); @test isnothing(e.calib)                  # Escape cancels
        ed.scale_bar!(e); ed.press!(e, 2.0, 1.0); ed.press!(e, 8.0, 1.0)                  # its ends: clicks, not edits
        @test length(e.calib) == 2 && !e.resizing
        ed.set_scale_length!(e, "-5"); ed.apply_scale!(e); @test b.width == 24.0            # no valid length: nothing
        ed.set_scale_length!(e, "5"); ed.apply_scale!(e)
        @test b.width ≈ 20.0 && ed.height(b) ≈ 10.0 && e.L ≈ [20.0, 10.0] && isnothing(e.calib)   # scaled, domain fitted
        ed.set_backdrop_field!(e, "x0", "-2"); @test b.x0 == -2.0 && e.L ≈ [20.0, 10.0]   # moved (cropped): domain kept
        ed.fit_domain!(e); @test e.L ≈ [18.0, 10.0]                                       # … until fitted again
        ed.set_backdrop_field!(e, "opacity", "3"); ed.set_backdrop_field!(e, "width", "abc")
        @test b.opacity == 1.0 && b.width ≈ 20.0                                            # clamped; ignored
        ed.remove_backdrop!(e); @test isnothing(e.backdrop) && occursin("\"backdrop\":null", ed.panel_json(e))
    end
    @testset "editing: snapping, shared edges, drag, undo" begin
        e = ed.Editor(Case(; L=[L1, L2], nel=[40, 10]))
        ed.new_zone!(e); foreach(p -> click!(e, p...), [(0.1, 0.1), (64.0, 0.05), (64.0, 4.0), (0.05, 4.05), (0.2, 0.2)])
        ed.new_zone!(e); foreach(p -> click!(e, p...), [(0.3, 4.3), (0.0, 12.0), (22.0, 12.0), (30.3, 4.3)]); ed.key!(e, "enter")
        A, B = e.regions[1].vertices, e.regions[2].vertices
        @test B[1] ∈ A && B[end][2] == A[3][2] && B[end] ∈ A                                # vertex snap; the edge snap split A's edge
        corner = copy(B[1])
        ed.move!(e, corner...); ed.press!(e, corner...); ed.move!(e, 0.0, 7*L2/20); ed.release!(e)
        @test e.regions[2].vertices[1] != corner && e.regions[2].vertices[1] ∈ e.regions[1].vertices   # dragged together
        ed.key!(e, "z"; cmd=true);  @test e.regions[2].vertices[1] == corner                # undo
        ed.key!(e, "z"; cmd=true, shift=true); @test e.regions[2].vertices[1] != corner     # redo
        ed.press!(e, e.regions[1].vertices[2]...; right=true)                              # right click: delete a vertex
        @test length(e.regions[1].vertices) == length(A) - 1
        e.ni = 3; ed.new_zone!(e); foreach(p -> click!(e, p...), [(40.1, 10.1), (50.0, 10.0), (45.0, 14.0)]); ed.key!(e, "enter")
        h3 = (L1/120, L2/30)
        @test all(v -> all(isapprox.(v ./ h3, round.(v ./ h3); atol=1e-6)), e.regions[3].vertices)
        ed.select_zone!(e, 3); ed.swap_zone!(e, -1); @test e.selected == 2
        ed.delete_zone!(e); @test length(e.regions) == 2
    end
    @testset "domain, zone constants, options, run, files" begin
        e = ed.Editor(Case(; L=[10.0, 5.0], nel=[20, 10]); step=1)
        ed.press!(e, 10.0, 5.0); ed.move!(e, 13.3, 6.7); ed.release!(e)
        @test e.L == [13.0, 7.0] && e.nel == [26, 14]                                       # snapped to 1 m, h = 0.5
        ed.set_domain!(e, e.L; h=1.0); @test e.nel == [13, 7]
        ed.set_bc!(e, "bottom", :fixed); @test e.bcs[2, 1] == :fixed
        ed.set_width!(e, "4"); ed.set_domain!(e, e.L); @test e.width == 4.0 && e.nely == 4
        ed.set_step!(e, 2); ed.new_zone!(e); foreach(p -> click!(e, p...), [(0, 0), (13, 0), (13, 3), (0, 3)]); ed.key!(e, "enter")
        ed.set_prop!(e, "E", "50"); ed.set_prop!(e, "nu", "0.2"); ed.set_prop!(e, "rho0", "2900")
        @test e.regions[1].props == Dict("E" => 5e7, "nu" => 0.2, "rho0" => 2900.0)
        which = findfirst(==("gimpm"), ed.option_entry("basis.which")[3])
        ed.set_option!(e, "basis.which", which); @test e.options["basis"]["which"] == "gimpm"
        ed.type_option!(e, "grf.param.σ", "7500"); @test e.options["grf"]["param"]["σ"] == 7500.0
        ed.type_option!(e, "grf.param.σ", "abc");  @test e.options["grf"]["param"]["σ"] == 7500.0   # invalid: ignored
        ed.tick_option!(e, "grf.param.σ", false);  @test !haskey(e.options, "grf")
        pv = ed.list_values("plot", ed.option_entry("plot.what")[3])
        ed.tick_value!(e, "plot.what", findfirst(v -> EP.option_label(v) == "plastic strain", pv), true)
        @test Set(EP.option_label.(e.options["plot"]["what"])) == Set(["pressure", "plastic strain"])
        ed.set_run!(e, "workflow", "elastoquasistatic!"); ed.set_run!(e, "lstps", "20"); ed.set_run!(e, "lstps", "-1")
        @test e.run["workflow"] == "elastoquasistatic!" && e.run["lstps"] == 20
        @test occursin("20 load steps", join(ed.checks(e).info))
        f = joinpath(mktempdir(), "q.toml"); ed.save!(e, f); @test isfile(f)
        pop!(e.regions); ed.load!(e, f); @test length(e.regions) == 1 && e.width == 4.0 && e.run["lstps"] == 20
        setoption!(d, s, p, v) = ed.setoption!(d, s, p, v)
        setoption!(e.options, "dtype", ["precision"], "32-bit precision (or single precision)")   # a warning: Done twice
        @test !ed.done!(e) && ed.done!(e)
        @test e.result.width == 4.0 && e.result.run.lstps == 20 && e.result.solver.basis.which == "gimpm"
        @test_throws ArgumentError ed.Editor(Case(; L=[1, 1], nel=[2, 2], zones=[Zone("r", Region(x -> x[1] < 0.5))]))
        @test ed.qml_index(nothing) == -1 && ed.julia_index(ed.qml_index(3)) == 3              # indices crossing to QML
        @test occursin("\"tab\":$(e.step - 1)", ed.panel_json(e))
        j = ed.json((; a=[1, 2.5], b="x\"y", c=nothing, d=true)); @test j == "{\"a\":[1.0,2.5],\"b\":\"x\\\"y\",\"c\":null,\"d\":true}"
    end
    @testset "cached points and checks, the cursor readout" begin
        e = ed.Editor(Case(; L=[L1, L2], nel=[40, 10]))
        ed.new_zone!(e); foreach(p -> click!(e, p...), [(0.0, 0.0), (64.0, 0.0), (64.0, 4.0), (0.0, 4.0)]); ed.key!(e, "enter")
        p = ed.preview_points(e)
        @test ed.preview_points(e) === p && p == ed.drawn_points(e.L, e.nel, e.ni, e.regions)  # cached, and right
        c = ed.checks(e); @test ed.checks(e) === c
        corner = copy(e.regions[1].vertices[3])
        ed.move!(e, corner...); ed.press!(e, corner...); ed.move!(e, 64.0, 8.0)
        @test ed.preview_points(e) === p                                                    # dragging: kept until release
        ed.release!(e)
        q = ed.preview_points(e)
        @test q !== p && q == ed.drawn_points(e.L, e.nel, e.ni, e.regions) && length(q[1]) > length(p[1])
        @test ed.checks(e) !== c                                                            # the zone changed
        c = ed.checks(e); ed.set_prop!(e, "E", "50"); @test ed.checks(e) !== c               # so did its constants
        status = e.status; ed.move!(e, 10.0, 2.0)
        @test e.status == status && startswith(e.readout, "x = ")                          # coordinates on the canvas
        @test occursin("\"readout\":\"x = ", ed.overlay_json(e, 800.0, 600.0))
        ed.leave!(e); @test e.readout == ""
    end

    @testset "without the window: how to get it" begin
        qt = EP.session_qt(); @test qt isa String                                           # Plots/GR load Qt
        for installed ∈ (false, true)
            m = EP.drawing_window_help(; installed)
            @test occursin("Qt6Base_jll@$qt", m) && occursin("_qt_version_tag", m)          # the pin, the symptom
            @test occursin("installed but not loaded", m) == installed
        end
        @test occursin("Qt6Base_jll@<the Qt version", EP.drawing_window_help(; installed=false, qt=nothing))
    end

    if isnothing(Base.find_package("QML"))
        @info "QML is not installed: the QML window itself is not tested"
    else
        ENV["QT_QPA_PLATFORM"] = get(ENV, "QT_QPA_PLATFORM", "offscreen"); ENV["QT_QUICK_BACKEND"] = "software"
        @eval using QML
        ext = Base.get_extension(EP, :DrawQMLExt)
        @testset "window defined" begin
            @test !isnothing(ext) && which(EP.drawing_window, Tuple{Case}).module == ext
        end
        @testset "the window, offscreen" begin
            d = mktempdir()
            case = Case(; L=[20.0, 10.0], nel=[20, 10], bcs=[:roller :roller; :fixed :free],
                         zones=[Zone("soil", Block([0, 0], [20, 6])), Zone("rock", Block([0, 0], [20, 2]); E=5e7)])
            @test isnothing(ext.qml_window(case; shots=(d, [1, 2, 3, 4])))                     # closes itself; no Done: nothing
            @test all(k -> isfile(joinpath(d, "step$k.png")), 1:4)
            d2 = mktempdir()                                                                   # a second window, same session
            @test isnothing(ext.qml_window(case; shots=(d2, [2]))) && isfile(joinpath(d2, "step2.png"))
        end
    end
end
