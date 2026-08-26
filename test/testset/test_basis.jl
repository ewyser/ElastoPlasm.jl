using StaticArrays
using LaTeXStrings
import ElastoPlasm: eval_basis

"""
    test_basis(n::Integer=4, m::Integer=4; which=["bsmpm","gimpm","smpm","mlsmpm"], npts::Integer=300, save::Bool=true)

Build a 2D `n×m`-element mesh and, for every basis kind in `which`, sweep a single particle
along the mesh's mid-height node row (`row = n÷2+1`, 1-indexed out of the `(n+1)×(m+1)` node
grid, tracking every node in that row, `cols = 1:n+1`) and compute two things at each sweep
position:

1. **Shape value/gradient** — for each physical node at that row, its shape value `Nᵢ(x)`
   and x-gradient `∂Nᵢ/∂x(x)` (`eval_basis`).
2. **Partition of unity** — `Σᵢ Nᵢ(x)`, summed over every candidate neighbor, which should
   equal 1 everywhere for a consistent basis.

Returns `(fig_shape, fig_pou)` — two `LaTeXStrings`-labelled figures, each with one subplot
per basis kind: `fig_shape` overlays the mid-height row's nodes' `Nᵢ` (solid) / `∂Nᵢ/∂x`
(dashed) curves; `fig_pou` plots `Σᵢ Nᵢ(x)` against the `Σᵢ Nᵢ = 1` reference line. Also
prints `max|ΣNᵢ-1|` per basis kind (informational, not a hard `@test` gate — some basis
kinds are not expected to hold exact PoU everywhere, e.g. this repo's uncorrected
`GimpBasis` kernel near boundaries; see `CLAUDE.md`'s Known bugs section).

# Arguments
- `n::Integer`, `m::Integer`: element counts along x/y (default 4×4).
- `which::Vector{String}`: basis kinds to test.
- `npts::Integer`: number of sweep positions along the line.
- `save::Bool`: if true, saves the figures under `dump/test/basis_shape/` and `basis_pou/`.

# Example
```julia
fig_shape, fig_pou = test_basis(4, 4)
```
"""
function test_basis(n::Integer=4, m::Integer=4; which::Vector{String}=["bsmpm","gimpm","smpm","mlsmpm"],
                     npts::Integer=300, save::Bool=true)
    cols   = 1:n+1
    row    = n÷2 + 1

    L, nel = Float64.([n, m]), [n, m]
    hx, hy = L[1]/n, L[2]/m
    target_x = [hx*(c-1) for c ∈ cols]
    ysweep   = hy*(row-1)
    xs       = range(1e-6, L[1]-1e-6, length=npts)

    shape_panels = Any[]
    pou_panels   = Any[]

    @testset "+ test_basis($n,$m)" verbose = true begin
    for w ∈ which
      @testset "- basis/$(w)" begin
        jld2  = slump_problem(L, nel; basis=(;get_default().basis...,which=w,how="Uii"), plot=(;status=false,freq=1.0,dpi=300,what=[]), fid="test/basis_$(w)")
        problem = load(jld2, "ic/problem")
        mesh  = problem.mesh
        mpts  = problem.mpts
        basis = load(jld2, "ic/basis")
        solver= load(jld2, "cfg/solver")

        nxnodes    = mesh.prprt.nno[1]
        target_ids = Int64[1 + (c-1) + nxnodes*(row-1) for c ∈ cols]
        NN = size(basis.p2n[1])[1]
        ip = Int64(1)

        Nline   = zeros(npts, length(cols))
        ∂Nline  = zeros(npts, length(cols))
        pouline = zeros(npts)
        for (k,x) ∈ enumerate(xs)
            mpts.x[ip] = SVector{2,Float64}(x, ysweep)
            solver.cairn.ignite.tplgy!(mpts, mesh, basis; ndrange=1); sync(CPU())
            s = 0.0
            for nn ∈ 1:NN
                no, N, ∂N = eval_basis(mpts, mesh, basis, ip, Int64(nn))
                s += N
                idx = findfirst(==(no), target_ids)
                if idx !== nothing
                    Nline[k,idx]  = N
                    ∂Nline[k,idx] = ∂N[1]
                end
            end
            pouline[k] = s
        end

        # 1) shape value/gradient
        pnl_shape = plot(xs, Nline; label=false, lw=1.5, ls=:solid,
                         title=w, titlefontsize=10, legend=false,
                         xlabel=L"x", ylabel=L"N_i(x),\ \partial N_i/\partial x",
                         ylims=(-1,1), framestyle=:box)
        plot!(pnl_shape, xs, ∂Nline; label=false, lw=1.0, ls=:dash)
        vline!(pnl_shape, target_x; lc=:gray, ls=:dot, lw=0.75, label=false)
        push!(shape_panels, pnl_shape)

        # 2) partition of unity, informational (not a hard @test gate)
        err    = maximum(abs.(pouline .- 1.0))
        logerr = log10(max(err, eps(Float64)))
        println("  $w: max|ΣN-1| along y=$(round(ysweep,digits=2)) = $err  (log10 = $logerr)")

        title_pou = LaTeXString("$(w):  \$\\log_{10}\\left(\\max|\\Sigma_i N_i - 1|\\right) = $(round(logerr,digits=2))\$")
        pnl_pou = plot(xs, pouline; label=false, lw=1.5,
                       title=title_pou, titlefontsize=10,
                       xlabel=L"x", ylabel=L"\sum_i N_i(x)",
                       ylims=(0.95,1.05), framestyle=:box)
        hline!(pnl_pou, [1.0]; lc=:gray, ls=:dash, lw=1.0, label=false)
        vline!(pnl_pou, target_x; lc=:gray, ls=:dot, lw=0.75, label=false)
        push!(pou_panels, pnl_pou)
      end
    end
    end

    ncols = length(which) > 1 ? 2 : 1
    nrows = ceil(Int, length(which)/ncols)

    # shared legend panel (one entry per tracked node, colors matching each subplot's
    # default color cycle) placed outside the grid instead of inside any one subplot
    legend_panel = plot(; legend=:left, grid=false, showaxis=false, ticks=nothing, framestyle=:none)
    for (i,c) ∈ enumerate(cols)
        plot!(legend_panel, [NaN], [NaN]; label=latexstring("N_{$row,$c}(x)"), lw=1.5, lc=i)
    end
    l = @layout [grid(nrows,ncols) a{0.16w}]
    fig_shape = plot(shape_panels..., legend_panel; layout=l, size=(560*ncols+140,380*nrows))
    fig_pou   = plot(pou_panels...;   layout=(nrows,ncols), size=(480*ncols,360*nrows))
    if save
        dir_shape = joinpath(ElastoPlasm.self.sys.dump, "test", "basis_shape")
        dir_pou   = joinpath(ElastoPlasm.self.sys.dump, "test", "basis_pou")
        mkpath(dir_shape); mkpath(dir_pou)
        savefig(fig_shape, joinpath(dir_shape, "basis_shape.png"))
        savefig(fig_pou,   joinpath(dir_pou,   "basis_pou.png"))
    end
    return fig_shape, fig_pou
end

test_basis(4, 4)
