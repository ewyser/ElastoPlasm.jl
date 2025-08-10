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
        #=
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
        =#
        return err
    end

    cases  = [
        (; which = "bsmpm", how = nothing     , ghost = false),
        #(; which = "gimpm", how = "undeformed", ghost = true ),
        #(; which = "smpm" , how = nothing     , ghost = true ),
    ]
    viz = (; status=false, freq=1.0, what=["P"], dims=(500.0,250.0) )

    nels = [[5, 10],[5, 20]]
    nk   = length(nels)+1;
    test = (; error = zeros(nk), h = zeros(nk), filename = [],)
    test.error[1],test.h[1] = Inf,Inf
    for basis ∈ cases
        @info "Testing with $(basis.which) basis"
        # 2d elastic collapse tests
        for (k,nel) ∈ enumerate(nels)
            @testset "- 2d geometry with $(basis.which) basis and nel = $nel" verbose = true begin
                l0      = 50.0
                ic, cfg = ic_collapse(nel, 0.0, 1.0e4, 80.0, l0; fid = "test/collapse", plot = viz, basis = basis,)
                err     = elastic_collapse(ic,cfg,l0)
                test.error[k+1] = err
                test.h[k+1]     = minimum(ic.mesh.h)
                test.filename
                @test test.error[k+1] < test.error[k] 

                if k == nk-1
                    @info "fdsfsdfdsfdsfsdfsdfsdfsd"
                    opts = (;
                        dims    = (250,250),
                        backend = gr(legend=true,markersize=2.0,markershape=:circle,markerstrokewidth=0.75,),
                        tit     = "Unnamed",
                        file    = joinpath(cfg[:paths][:plot],"2d_elastic_column_collapse_convergence.png"),
                    )
                    config_plot()
                    opts.backend

                    p = plot(
                        test.h,test.error,
                        seriestype =:line, 
                        xlabel     = L"$x-$direction"*" [m]",
                        ylabel     = L"$z-$direction"*" [m]",
                        label      ="$(length(nel))D $(cfg[:instr][:basis][:which]), $(cfg[:instr][:fwrk][:trsfr]) mapping",
                        title      = "$(opts.tit)",
                        size       = opts.dims,
                    )
                    display(plot(p; layout=(1,1), size=(450,250)))
                    save_plot(opts)
                end
            end
        end

        #==#
    end


        # 3d slump tests
        # TODO: Add 3d elastic collapse tests


end