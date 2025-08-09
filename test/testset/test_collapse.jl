@testset "+ $(basename(@__FILE__))" verbose = true begin
    function elastic_collapse(ic,cfg,l0)
        dim     = size(ic.mpts.x,1)
        z0      = copy(ic.mpts.x[end, :])
        out     = collapse(ic, cfg)
        mesh,mpts,cmpr = out.ic.mesh, out.ic.mpts, out.ic.cmpr

        # Numeric and analytic solution
        idx = mesh.dim == 2 ? 2 : 3
        xnum = abs.(mpts.s.σᵢ[idx, :])
        ynum = z0
        x = abs.(cmpr.ρ0 * 9.81 * (l0 .- z0))
        y = z0
        err = sum(sqrt.((xnum .- x).^2) .* mpts.Ω₀) / (9.81 * cmpr.ρ0 * l0 * sum(mpts.Ω₀))
        
        config_plot()
        gr(size=(500,300),legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.0,markerstrokecolor=:match,)
        p1 = plot(
            xnum .* 1e-3, ynum,
            seriestype = :scatter,
            label = "$(dim)d $(cfg.instr[:basis][:which])",
            markersize = 2.0,
            xlabel = L"\sigma_{yy}\ \mathrm{[kPa]}",
            ylabel = L"y\ \mathrm{[m]}",
            legend = :topright,
            grid = true,
            color = :red,
            dpi = 120
        )
        plot!(
            x .* 1e-3, y,
            seriestype = :line,
            label = L"\sum_{p}\frac{||\sigma_{yy}^p - \sigma_{yy}^a(x_p)|| \Omega_0^p}{(g \rho_0 l_0) \Omega_0}", 
            linewidth = 2,
            color = :blue
        )
        display(plot(p1; layout = (1, 1), size = (500, 300), dpi = 120))
        savefig(joinpath(cfg.paths[:plot],"$(dim)d_numeric_analytic_$(basename(@__FILE__)).png"))
        return err
    end

    cases  = [
        (; which = "bsmpm", how = nothing     , ghost = false),
        #(; which = "gimpm", how = "undeformed", ghost = true ),
        #(; which = "smpm" , how = nothing     , ghost = true ),
    ]
    viz = (; status=true, freq=1.0, what=["P"], dims=(500.0,250.0) )

    Error = zeros(4); Error[1] = Inf

    for basis ∈ cases
        @info "Testing with $(basis.which) basis"
        # 2d elastic collapse tests
        for (k,nel) ∈ enumerate([[5, 10],[5, 20],[5, 40]])
            @testset "- 2d geometry with $(basis.which) basis and nel = $nel" verbose = true begin
                l0      = 50.0
                ic, cfg = ic_collapse(nel, 0.0, 1.0e4, 80.0, l0; fid = "test/collapse", plot = viz, basis = basis,)
                Error[k+1] = elastic_collapse(ic,cfg,l0)
                @test Error[k+1] < Error[k] 
            end
        end
        

        # 3d slump tests
        # TODO: Add 3d elastic collapse tests
    end





end