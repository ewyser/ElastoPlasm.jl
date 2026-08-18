using StaticArrays
import ElastoPlasm: eval_basis

"""
    test_basis(n::Integer=4, m::Integer=4; which=["bsmpm","gimpm","smpm","mlsmpm"], npts::Integer=300, nx::Integer=60, ny::Integer=60, row::Integer=4, cols=1:4, save::Bool=true)

Build a 2D `n×m`-element mesh and, for every basis kind in `which`, compute two things:

1. **Shape value/gradient along a line** — sweep a particle along `y = (row-1)*h[2]` (the
   physical node row selected by `row`) and track, for each of the physical nodes at that
   row/`cols` (1-indexed, out of the `(n+1)×(m+1)` node grid), its shape value `N(x)` and
   x-gradient `∂N/∂x(x)`.
2. **Partition of unity across the whole mesh** — sweep a particle over a fine 2D grid
   covering the *entire* mesh and compute `ΣNᵢ(xp)` at every point (mirroring the
   `dev`-branch mlsmpm prototype's 2D heatmap verification, generalized from "one node's
   kernel" to "PoU everywhere").

Returns `(fig_shape, fig_pou)` — two figures, each with one subplot per basis kind (not per
node): `fig_shape`'s panels overlay the `row`/`cols` nodes' `N`/`∂N/∂x` curves; `fig_pou`'s
panels are PoU heatmaps with those same nodes marked as red crosses for reference. Also
prints `max|ΣN-1|` per basis kind (informational, not a hard `@test` gate — some basis kinds
are not expected to hold exact PoU everywhere, e.g. this repo's uncorrected `GimpBasis`
kernel near boundaries; see `CLAUDE.md`'s Known bugs section).

# Arguments
- `n::Integer`, `m::Integer`: element counts along x/y (default 4×4).
- `which::Vector{String}`: basis kinds to test.
- `npts::Integer`: number of sweep positions for the line sweep (1).
- `nx::Integer`, `ny::Integer`: sweep-grid resolution along x/y for the PoU heatmap (2).
- `row::Integer`, `cols`: 1-indexed physical node-grid coordinates to track/mark.
- `save::Bool`: if true, saves the figures under `dump/test/basis_shape/` and `basis_pou/`.

# Example
```julia
fig_shape, fig_pou = test_basis(4, 4)
```
"""
function test_basis(n::Integer=4, m::Integer=4; which::Vector{String}=["bsmpm","gimpm","smpm","mlsmpm"],
                     npts::Integer=300, nx::Integer=60, ny::Integer=60, row::Integer=4, cols=1:4, save::Bool=true)
    L, nel = Float64.([n, m]), [n, m]
    hx, hy = L[1]/n, L[2]/m
    target_x  = [hx*(c-1) for c ∈ cols]
    target_y  = fill(hy*(row-1), length(cols))
    ysweep    = hy*(row-1)

    xs_line = range(1e-6, L[1]-1e-6, length=npts)
    xs_grid = range(1e-6, L[1]-1e-6, length=nx)
    ys_grid = range(1e-6, L[2]-1e-6, length=ny)

    shape_panels = Any[]
    pou_panels   = Any[]

    @testset "+ test_basis($n,$m)" verbose = true begin
    for w ∈ which
      @testset "- basis/$(w)" begin
        jld2  = ic_slump(L, nel; basis=(;which=w,how="Uii",ghost=0), plot=(;status=false,freq=1.0,dpi=300,what=[]), fid="test/basis_$(w)")
        mesh  = load(jld2, "ic/mesh")
        mpts  = load(jld2, "ic/mpts")
        basis = load(jld2, "ic/basis")
        solver= load(jld2, "cfg/solver")

        nxnodes    = mesh.prprt.nno[1]
        target_ids = Int64[1 + (c-1) + nxnodes*(row-1) for c ∈ cols]
        NN = size(basis.p2n[1])[1]
        ip = Int64(1)

        # 1) shape value/gradient of the target nodes along the sweep line
        Nline  = zeros(npts, length(cols))
        ∂Nline = zeros(npts, length(cols))
        for (k,x) ∈ enumerate(xs_line)
            mpts.x[ip] = SVector{2,Float64}(x, ysweep)
            solver.cairn.ignite.tplgy!(mpts, mesh, basis; ndrange=1); sync(CPU())
            for nn ∈ 1:NN
                no, N, ∂N = eval_basis(mpts, mesh, basis, ip, Int64(nn))
                idx = findfirst(==(no), target_ids)
                if idx !== nothing
                    Nline[k,idx]  = N
                    ∂Nline[k,idx] = ∂N[1]
                end
            end
        end
        labels = reshape(["node[$(row),$(c)]" for c ∈ cols], 1, :)
        pnl_shape = plot(xs_line, Nline; label=labels, lw=1.5, ls=:solid,
                         title="$w", titlefontsize=9, legend=(w==first(which) ? :outerright : false))
        plot!(pnl_shape, xs_line, ∂Nline; label=false, lw=1.0, ls=:dash)
        push!(shape_panels, pnl_shape)

        # 2) partition of unity over the whole mesh
        pou = zeros(ny, nx)
        for (j,y) ∈ enumerate(ys_grid), (i,x) ∈ enumerate(xs_grid)
            mpts.x[ip] = SVector{2,Float64}(x, y)
            solver.cairn.ignite.tplgy!(mpts, mesh, basis; ndrange=1); sync(CPU())
            s = 0.0
            for nn ∈ 1:NN
                _, N, _ = eval_basis(mpts, mesh, basis, ip, Int64(nn))
                s += N
            end
            pou[j,i] = s
        end
        err = maximum(abs.(pou .- 1.0))
        println("  $w: max|ΣN-1| over full $(n)x$(m)-element mesh = $err")

        pnl_pou = heatmap(xs_grid, ys_grid, pou; title="$w (max|ΣN-1|=$(round(err,sigdigits=3)))", titlefontsize=9,
                          c=:viridis, clims=(0.9,1.1), xlabel="x", ylabel="y", aspect_ratio=:equal)
        scatter!(pnl_pou, target_x, target_y; m=:xcross, ms=6, mc=:red, label=false)
        push!(pou_panels, pnl_pou)
      end
    end
    end

    ncols = length(which) > 1 ? 2 : 1
    nrows = ceil(Int, length(which)/ncols)
    fig_shape = plot(shape_panels...; layout=(nrows,ncols), size=(560*ncols,380*nrows),
                     plot_title="N (solid) / ∂N/∂x (dash) along y=$(round(ysweep,digits=2)), nodes[row=$row,cols=$cols]",
                     plot_titlefontsize=10)
    fig_pou = plot(pou_panels...; layout=(nrows,ncols), size=(520*ncols,440*nrows),
                   plot_title="Partition of unity across $(n)x$(m)-element mesh", plot_titlefontsize=11)
    if save
        dir_shape = joinpath(ElastoPlasm.self.sys.out, "test", "basis_shape")
        dir_pou   = joinpath(ElastoPlasm.self.sys.out, "test", "basis_pou")
        mkpath(dir_shape); mkpath(dir_pou)
        savefig(fig_shape, joinpath(dir_shape, "basis_shape.png"))
        savefig(fig_pou,   joinpath(dir_pou,   "basis_pou.png"))
    end
    return fig_shape, fig_pou
end

test_basis(4, 4)
