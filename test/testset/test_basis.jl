using StaticArrays
import ElastoPlasm: eval_basis

"""
    test_basis(n::Integer, m::Integer; which=["bsmpm","gimpm","smpm","mlsmpm"], npts::Integer=300, save::Bool=true)

Build a 2D `n×m`-element mesh and, for every basis kind in `which`, sweep a single material
point along the mesh's mid-height horizontal line and evaluate every candidate neighbor
node's shape value/gradient (`eval_basis`) at each sweep position.

At each sweep position, checks:
- partition of unity:    `Σ Nᵢ = 1`
- linear consistency:    `Σ Nᵢ·δxᵢ = 0`   (δxᵢ = node position - particle position)
- gradient consistency:  `Σ ∂Nᵢ = 0`

reported (printed, and annotated on each figure's title) rather than hard-gated with `@test`
— some basis kinds are not expected to hold these exactly everywhere (this repo's uncorrected
`GimpBasis` kernel does not reproduce exact PoU near boundaries, a known, pre-existing
property; see `CLAUDE.md`'s Known bugs section). Returns `Dict{String,Plots.Plot}` — one
figure per basis kind, one panel per candidate neighbor-node index, showing that node's
shape value `N(xp)` and x-gradient `∂N/∂x(xp)` swept over the sweep line. Node-panel plots
double as a visual check too (e.g. `BSplineBasis`'s piecewise boundary correction should show
up as kinks at the node-type transition, `MLSBasis` should stay smooth everywhere including
boundaries).

# Arguments
- `n::Integer`, `m::Integer`: element counts along x/y.
- `which::Vector{String}`: basis kinds to test (`get_option().basis.which` lists valid values).
- `npts::Integer`: number of sweep positions along the mid-height line.
- `save::Bool`: if true, also saves each figure to `dump/test/basis_<kind>/basis_<kind>.png`.

# Example
```julia
figs = test_basis(6, 4)
figs["mlsmpm"]  # the MLS figure
```
"""
function test_basis(n::Integer, m::Integer; which::Vector{String}=["bsmpm","gimpm","smpm","mlsmpm"], npts::Integer=300, save::Bool=true)
    L, nel = Float64.([n, m]), [n, m]
    figs = Dict{String,Any}()

    @testset "+ test_basis($n,$m)" verbose = true begin
    for w ∈ which
      @testset "- basis/$(w)" begin
        jld2  = ic_slump(L, nel; basis=(;which=w,how="Uii",ghost=0), plot=(;status=false,freq=1.0,dpi=300,what=[]), fid="test/basis_$(w)")
        mesh  = load(jld2, "ic/mesh")
        mpts  = load(jld2, "ic/mpts")
        basis = load(jld2, "ic/basis")
        solver= load(jld2, "cfg/solver")

        NN   = size(basis.p2n[1])[1]
        ymid = mesh.prprt.L[2] / 2
        xs   = collect(range(1e-6, mesh.prprt.L[1] - 1e-6, length=npts))

        N   = zeros(npts, NN)
        ∂Nx = zeros(npts, NN)
        ∂Ny = zeros(npts, NN)
        pou    = zeros(npts)
        lincon = zeros(npts, 2)
        gradcon= zeros(npts, 2)

        ip = Int64(1)
        for (k,x) ∈ enumerate(xs)
            mpts.x[ip] = SVector{2,Float64}(x, ymid)
            solver.cairn.ignite.tplgy!(mpts, mesh, basis; ndrange=1); sync(CPU())

            active = Int64[]
            for nn ∈ 1:NN
                no, Nv, ∂Nv = eval_basis(mpts, mesh, basis, ip, Int64(nn))
                N[k,nn]   = Nv
                ∂Nx[k,nn] = ∂Nv[1]
                ∂Ny[k,nn] = ∂Nv[2]
                if !iszero(no)
                    push!(active, nn)
                end
            end
            pou[k] = sum(N[k,:])
            δx = [mesh.x[basis.p2n[ip][nn]] .- mpts.x[ip] for nn ∈ active]
            lincon[k,:]  .= sum(N[k,active[i]].*δx[i] for i ∈ eachindex(δx))
            gradcon[k,:] .= (sum(∂Nx[k,:]), sum(∂Ny[k,:]))
        end

        # informational, not a hard pass/fail gate (matches test_performance.jl's style) —
        # some basis kinds are not expected to hold these exactly everywhere (e.g. this
        # repo's uncorrected GimpBasis kernel does not reproduce exact PoU near boundaries;
        # see CLAUDE.md's known-bugs section for basis-kind-specific caveats)
        pou_err  = maximum(abs.(pou .- 1.0))
        lin_err  = maximum(abs.(lincon))
        grad_err = maximum(abs.(gradcon))
        println("  $w: max|ΣN-1|=$pou_err  max|Σ N·δx|=$lin_err  max|Σ∂N|=$grad_err")

        ncols = ceil(Int, sqrt(NN))
        nrows = ceil(Int, NN/ncols)
        panels = [
            plot(xs, [N[:,nn] ∂Nx[:,nn]];
                 label=["N" "∂N/∂x"], title="node $nn", titlefontsize=8,
                 legend = nn==1 ? :topright : false, lw=1.5)
            for nn ∈ 1:NN
        ]
        title = "$w — shape value/gradient vs particle x (y=$(round(ymid,digits=2)))\n" *
                "max|ΣN-1|=$(round(pou_err,sigdigits=3))  max|ΣN·δx|=$(round(lin_err,sigdigits=3))  max|Σ∂N|=$(round(grad_err,sigdigits=3))"
        fig = plot(panels...; layout=(nrows,ncols), size=(300*ncols,220*nrows+40),
                   plot_title=title, plot_titlefontsize=10)
        figs[w] = fig
        if save
            dir = joinpath(ElastoPlasm.self.sys.out, "test", "basis_$(w)")
            mkpath(dir)
            savefig(fig, joinpath(dir, "basis_$(w).png"))
        end
      end
    end
    end
    return figs
end

test_basis(6, 4)
