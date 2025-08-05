@testset "+ $(basename(@__FILE__))" verbose = true begin

    nel  = [5, 10]
    dim  = length(nel)
    ν, E, ρ0, l0 = 0.0, 1.0e4, 80.0, 50.0
    ic, cfg = ic_collapse(nel, ν, E, ρ0, l0; fid = "test/collapse", plot = (; status=true, freq=1.0, what=["P"], dims=(500.0,250.0) ))
    z0      = copy(ic.mp.x[end, :])
    out     = collapse(ic, cfg)
    mesh,mp,cmpr = out.ic.mesh, out.ic.mp, out.ic.cmpr

    # Numeric and analytic solution
    idx = mesh.dim == 2 ? 2 : 3
    xnum = abs.(mp.s.σᵢ[idx, :])
    ynum = z0
    x = abs.(cmpr.ρ0 * 9.81 * (l0 .- z0))
    y = z0

    err = sum(sqrt.((xnum .- x).^2) .* mp.Ω₀) / (9.81 * cmpr.ρ0 * l0 * sum(mp.Ω₀))
    
    @test isa(err,AbstractFloat)

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
end