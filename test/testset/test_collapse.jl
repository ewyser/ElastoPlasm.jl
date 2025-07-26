@testset "+ $(basename(@__FILE__))" verbose = true begin

    nel  = [5, 10]
    dim  = length(nel)
    ν, E, ρ0, l0 = 0.0, 1.0e4, 80.0, 10.0
    ic, cfg = ic_collapse(nel, ν, E, ρ0, l0; fid = "test/collapse", plot = (; status=true, freq=1.0, what=["P"], dims=(500.0,250.0) ))
    z0      = copy(ic.mp.x[end, :])
    out     = collapse(ic, cfg)
    mesh,mp,cmpr = out.mesh, out.mp, ic.cmpr

    # Numeric and analytic solution
    idx = mesh.dim == 2 ? 2 : 3
    xnum = abs.(mp.σᵢ[idx, :])
    ynum = z0
    x = abs.(cmpr.ρ0 * 9.81 * (l0 .- z0))
    y = z0

    err = sum(sqrt.((xnum .- x).^2) .* mp.Ω₀) / (9.81 * cmpr.ρ0 * l0 * sum(mp.Ω₀))
    #println(err)

    config_plot()
    #=    =#
    p1 = scatter(xnum .* 1e-3, ynum, label="$(dim)")
    #plot!(x .* 1e-3, y, label=L"\sum_{p}\dfrac{||\sigma_{yy}^p-\sigma_{yy}^a(x_p)||V_0^p}{(g\rho_0l_0)V_0}", xlabel=L"$\sigma_{yy}$ [kPa]", ylabel=L"$y-$position [m]")
    display(plot(p1; layout=(1,1), size=(450,250)))
    savefig(joinpath(cfg.paths[:plot],"$(dim)d_numeric_analytic_$(basename(@__FILE__)).png"))
end