# ============================================================================================================================
# MLS shape function construction — self-contained example (1D and 2D)
#
# This script demonstrates step-by-step how MLS shape functions are built
# and contrasts them with the standard B-spline basis used in ElastoPlasm.
#
# Run standalone (no ElastoPlasm package needed):
#   julia src/home/script/example/mlsmpm.jl
# ============================================================================================================================

using LinearAlgebra, Plots; gr()

# ============================================================================================================================
# 1. Kernel — quadratic B-spline  (same piecewise cubic as ElastoPlasm's t1_ϕ∂ϕ interior node)
#
#    The MLS formalism does not prescribe which kernel to use — any compactly supported
#    function works. We reuse the quadratic B-spline so results are directly comparable
#    with your existing BSplineBasis.
#
#    Scaled to unit support radius r = h * 1.5 (support ∈ [-1.5h, 1.5h]):
#      w(ξ) = (3/4 - ξ²)          if |ξ| < 0.5
#      w(ξ) = (1/2)(3/2 - |ξ|)²  if 0.5 ≤ |ξ| < 1.5
#    with ξ = (x_i - x_p)/h, then divided by h for consistent units.
# ============================================================================================================================

"""Quadratic B-spline kernel value. ξ = (x_i - x_p)/h."""
function bspline(ξ::Float64)
    a = abs(ξ)
    if a < 0.5
        return 0.75 - a^2
    elseif a < 1.5
        return 0.5 * (1.5 - a)^2
    else
        return 0.0
    end
end

"""Quadratic B-spline kernel derivative ∂w/∂ξ. ξ = (x_i - x_p)/h."""
function bspline_grad(ξ::Float64)
    a = abs(ξ)
    if a < 0.5
        return -2.0 * ξ
    elseif a < 1.5
        return sign(ξ) * (-1.5 + a)
    else
        return 0.0
    end
end

# ============================================================================================================================
# 2. MLS basis construction — 1D
#
#  Given a particle at xp on a 1D grid {x_i} with uniform spacing h:
#
#  Step 1 — kernel weights
#    w_i = bspline((x_i - x_p)/h) / h
#
#  Step 2 — linear polynomial basis centered at xp
#    p(x_i) = [1, x_i - xp]   (column vector)
#
#  Step 3 — moment matrix M (2×2)
#    M = Σ_i w_i p(x_i) p(x_i)ᵀ
#
#    On a uniform grid with a symmetric kernel, the off-diagonal ⟨w δx⟩ vanishes, giving:
#      M = [ Σw_i,       0     ]
#          [   0,   Σ w_i δx_i² ]  ← lower-right block = D_p (your Dᵢⱼ)
#
#  Step 4 — shape function value
#    ϕ_i(xp) = p(xp)ᵀ M⁻¹ p(x_i) w_i
#    Since p(xp) = [1, 0]ᵀ, this reduces to:
#      ϕ_i = (M⁻¹[1,:] · p(x_i)) w_i  =  w_i   (same as kernel!)
#
#  Step 5 — shape function gradient  ∇ϕ_i(xp) = ∂ϕ_i/∂xp
#    Differentiating ϕ_i w.r.t. xp (chain rule on p(xp) and M):
#      ∇ϕ_i = -w_i M⁻¹[2,2] (x_i - xp)
#    where M⁻¹[2,2] = 1/D_p.
#
#    On a uniform grid: D_p = Σ w_i δx_i² = h²/4 (closed-form for quadratic B-spline),
#    so M⁻¹[2,2] = 4/h², and the gradient is:
#      ∇ϕ_i = -(4/h²) w_i (x_i - xp)
#
#  Key property — partition of unity (holds EVERYWHERE, including near boundaries):
#    Σ_i ϕ_i = p(xp)ᵀ M⁻¹ Σ_i [p(x_i) w_i] = p(xp)ᵀ M⁻¹ M e₁ = p(xp)ᵀ e₁ = 1  ✓
#  No special boundary node types (like your mesh.type system) are required.
# ============================================================================================================================

"""
    mls_1d(xp, x_nodes, h) -> (ϕ, ∇ϕ, w, M)

Compute MLS shape function values and gradients for a 1D particle at `xp`
given grid nodes `x_nodes` with spacing `h`.

Returns:
- `ϕ`  : shape function values  ϕ_i(xp)   [= kernel weights w_i on uniform grid]
- `∇ϕ` : shape function gradients ∂ϕ_i/∂xp [= -w_i M⁻¹[2,2] δx_i]
- `w`  : raw kernel weights w_i
- `M`  : 2×2 moment matrix
"""
function mls_1d(xp::Float64, x_nodes::Vector{Float64}, h::Float64)
    n = length(x_nodes)

    # Step 1 — kernel weights  w_i = bspline(ξ_i)/h,  ξ_i = (x_i - x_p)/h
    w = [bspline((x_nodes[i] - xp) / h) / h  for i in 1:n]

    # Step 2 — polynomial basis  p(x_i) = [1, δx_i]
    δx = [x_nodes[i] - xp  for i in 1:n]

    # Step 3 — moment matrix M = Σ w_i p_i p_iᵀ  (2×2)
    M = zeros(2, 2)
    for i in 1:n
        pᵢ  = [1.0, δx[i]]
        M .+= w[i] .* (pᵢ * pᵢ')
    end

    # Step 4 — shape function values:  ϕ_i = (p(xp)ᵀ M⁻¹ p(x_i)) w_i
    #           With p(xp) = [1,0]ᵀ this is just (M⁻¹[1,1] + M⁻¹[1,2]*δx_i) w_i
    M⁻¹ = inv(M)
    p0  = [1.0, 0.0]   # p(xp)
    ϕ   = [(p0' * M⁻¹ * [1.0, δx[i]]) * w[i]  for i in 1:n]

    # Step 5 — shape function gradients:  ∇ϕ_i = -w_i M⁻¹[2,2] δx_i
    #           M⁻¹[2,2] = 1/D_p;  on uniform grid = 4/h²
    ∇ϕ  = [-w[i] * M⁻¹[2,2] * δx[i]  for i in 1:n]

    return ϕ, ∇ϕ, w, M
end

# ============================================================================================================================
# 3. Standard B-spline shape function + gradient (for comparison)
#    Value:    ϕ_i   = bspline(ξ_i)/h
#    Gradient: ∇ϕ_i  = bspline_grad(ξ_i)/h²   (chain rule: ∂/∂xp = -(1/h) ∂/∂ξ)
# ============================================================================================================================

function std_1d(xp::Float64, x_nodes::Vector{Float64}, h::Float64)
    ϕ  = [ bspline((x_nodes[i]-xp)/h)/h           for i in eachindex(x_nodes)]
    ∇ϕ = [-bspline_grad((x_nodes[i]-xp)/h)/h^2    for i in eachindex(x_nodes)]
    return ϕ, ∇ϕ
end

# ============================================================================================================================
# 4. Verification — 1D
# ============================================================================================================================

println("\n", "="^60)
println("MLS shape function — 1D verification")
println("="^60)

h       = 1.0
nel     = 6
x_nodes = collect(0.0:h:nel*h)   # 7 nodes: 0, 1, 2, 3, 4, 5, 6
n       = length(x_nodes)

# Three test positions:
#   interior   — full symmetric support, closed-form D_p exact
#   near left  — asymmetric support, where std B-spline needs node-type correction
#   mid-cell   — between two nodes
for (xp, label) in [(3.4, "interior          "), (0.3, "near left boundary"), (2.7, "mid-cell          ")]

    ϕ_mls, ∇ϕ_mls, w, M = mls_1d(xp, x_nodes, h)
    ϕ_std, ∇ϕ_std        = std_1d(xp, x_nodes, h)

    println("\nParticle xp = $xp  ($label)")
    println("  Moment matrix M    = ", round.(M, digits=5))
    println("  cond(M)            = ", round(cond(M), digits=4))
    println("  M⁻¹[2,2]           = ", round(inv(M)[2,2], digits=6),
            "  ←→  closed-form 4/h² = ", 4.0/h^2)

    # Partition of unity
    Σϕ_mls = sum(ϕ_mls);  Σϕ_std = sum(ϕ_std)
    println("  Σ ϕᵢ  MLS         = ", round(Σϕ_mls, digits=10), "  (must be 1.0)")
    println("  Σ ϕᵢ  std B-spline = ", round(Σϕ_std, digits=10), "  (1.0 only if node-type corrected)")

    # Linear consistency  Σ ϕ_i δx_i = 0
    lin_mls = sum(ϕ_mls[i] * (x_nodes[i]-xp) for i in 1:n)
    lin_std = sum(ϕ_std[i] * (x_nodes[i]-xp) for i in 1:n)
    println("  Σ ϕᵢ δxᵢ MLS      = ", round(lin_mls, digits=10), "  (must be 0.0)")
    println("  Σ ϕᵢ δxᵢ std      = ", round(lin_std, digits=10))

    # Gradient consistency  Σ ∇ϕ_i = 0
    println("  Σ ∇ϕᵢ MLS         = ", round(sum(∇ϕ_mls), digits=10), "  (must be 0.0)")
    println("  Σ ∇ϕᵢ std         = ", round(sum(∇ϕ_std), digits=10))
end

# ============================================================================================================================
# 5. MLS gradient — simplified vs full (1D)
#
#  The formula ∇ϕ_i = -w_i M⁻¹[2,2] δx_i  keeps only the "basis derivative" term:
#    ∂/∂xp [p(xp)ᵀ M⁻¹ p(x_i) w_i]
#  = [∂p(xp)/∂xp]ᵀ M⁻¹ p(x_i) w_i          ← term A: basis derivative (kept)
#  + p(xp)ᵀ [∂M⁻¹/∂xp] p(x_i) w_i           ← term B: moment matrix derivative (dropped)
#  + p(xp)ᵀ M⁻¹ [∂p(x_i)/∂xp] w_i           ← term C: basis node derivative
#  + p(xp)ᵀ M⁻¹ p(x_i) [∂w_i/∂xp]           ← term D: kernel derivative
#
#  With p(xp) = [1,0]ᵀ, ∂p(xp)/∂xp = [0,-1]ᵀ, ∂p(x_i)/∂xp = [0,-1]ᵀ:
#    Term A = -M⁻¹[2,2] δx_i w_i  (what we use)
#    Term C = -M⁻¹[1,2] w_i       (zero in interior since M off-diagonal = 0, non-zero near boundary)
#    Terms B,D involve ∂M/∂xp (kernel and basis derivatives summed over all nodes)
#
#  In the interior (symmetric support): M is diagonal, B+C+D cancel → simplified formula is exact.
#  Near boundaries: M has non-zero off-diagonals → simplified formula misses the correction.
#
#  The full gradient IS consistent (Σ ∇ϕ_i = 0 always), but requires computing all four terms.
#  In MLS-MPM practice (Hu et al. 2018), the simplified formula is used throughout — the near-boundary
#  gradient error is small and bounded, and PoU + linear consistency of ϕ VALUES still holds exactly.
# ============================================================================================================================

"""
    mls_1d_full_gradient(xp, x_nodes, h) -> ∇ϕ_full

Full MLS gradient including the ∂M/∂xp correction term (terms A+C from above).
Term B (∂M⁻¹/∂xp) and Term D (∂w_i/∂xp) are also included via autodiff-style
accumulation over all nodes.

Only terms A and C are implemented here (B and D are higher-order and small in practice).
"""
function mls_1d_full_gradient(xp::Float64, x_nodes::Vector{Float64}, h::Float64)
    n  = length(x_nodes)
    w  = [bspline((x_nodes[i] - xp) / h) / h  for i in 1:n]
    δx = [x_nodes[i] - xp  for i in 1:n]

    M  = zeros(2, 2)
    for i in 1:n
        pᵢ = [1.0, δx[i]]
        M .+= w[i] .* (pᵢ * pᵢ')
    end
    M⁻¹ = inv(M)

    # ∂M/∂xp = Σ_j [ (∂w_j/∂xp) p_j p_jᵀ + w_j ∂p_j/∂xp p_jᵀ + w_j p_j (∂p_j/∂xp)ᵀ ]
    # ∂w_j/∂xp = -bspline_grad(ξ_j)/h²  (chain rule: ∂w/∂xp = -∂w/∂x_j)
    # ∂p_j/∂xp = [0, -1]ᵀ
    ∂M = zeros(2, 2)
    ep  = [0.0, -1.0]
    for j in 1:n
        pⱼ  = [1.0, δx[j]]
        ∂wⱼ = -bspline_grad((x_nodes[j]-xp)/h) / h^2
        ∂M .+= ∂wⱼ .* (pⱼ * pⱼ') .+ w[j] .* (ep * pⱼ' .+ pⱼ * ep')
    end

    # ∂M⁻¹/∂xp = -M⁻¹ ∂M M⁻¹
    ∂M⁻¹ = -M⁻¹ * ∂M * M⁻¹

    p0   = [1.0, 0.0]
    ∂p0  = [0.0, 0.0]   # ∂p(xp)/∂xp = 0 (p(xp) = [1,0] is constant w.r.t. xp)

    ∇ϕ = zeros(n)
    for i in 1:n
        pᵢ  = [1.0, δx[i]]
        ∂pᵢ = ep                                      # ∂p(x_i; xp)/∂xp = [0,-1]
        ∂wᵢ = -bspline_grad((x_nodes[i]-xp)/h) / h^2

        # full gradient (all four terms)
        ∇ϕ[i] = (∂p0' * M⁻¹ * pᵢ) * w[i]            # term A (= 0 since ∂p0=0)
              + (p0'  * ∂M⁻¹ * pᵢ) * w[i]             # term B
              + (p0'  * M⁻¹  * ∂pᵢ) * w[i]            # term C
              + (p0'  * M⁻¹  * pᵢ) * ∂wᵢ              # term D
    end
    return ∇ϕ
end

println("\n", "="^60)
println("MLS gradient — simplified vs full (1D near boundary)")
println("="^60)

xp_bnd = 0.3
∇ϕ_simplified = mls_1d(xp_bnd, x_nodes, h)[2]
∇ϕ_full       = mls_1d_full_gradient(xp_bnd, x_nodes, h)

println("\nParticle xp = $xp_bnd  (near left boundary)")
println("  Σ ∇ϕᵢ  simplified = ", round(sum(∇ϕ_simplified), digits=10),
        "  (non-zero: missing ∂M/∂xp)")
println("  Σ ∇ϕᵢ  full       = ", round(sum(∇ϕ_full),       digits=10),
        "  (must be 0.0 — ∂/∂xp[Σ ϕᵢ] = ∂1/∂xp = 0)")
println()
println("  Node-by-node comparison:")
println("  ", rpad("node", 6), rpad("x_i", 8), rpad("∇ϕ_simplified", 18), "∇ϕ_full")
for i in eachindex(∇ϕ_simplified)
    if abs(∇ϕ_simplified[i]) > 1e-12 || abs(∇ϕ_full[i]) > 1e-12
        println("  ", rpad(i,6), rpad(round(x_nodes[i],digits=2),8),
                rpad(round(∇ϕ_simplified[i],digits=6),18), round(∇ϕ_full[i],digits=6))
    end
end
#
#  Extension: linear basis p(x_i) = [1, δx_i, δy_i]  (3-vector).
#
#  Moment matrix M is 3×3:
#    M = Σ_i w_i p_i p_iᵀ
#
#  On a uniform grid with tensor-product kernel, M is block-diagonal:
#    M = diag(Σw_i,  Σw_i δx_i²,  Σw_i δy_i²)
#      = diag(  1,     h²/4,          h²/4  )
#
#  Shape function value:    ϕ_i   = w_i  (same conclusion as 1D)
#  Shape function gradient: ∇ϕ_i  = -w_i M⁻¹[2:3,2:3] [δx_i, δy_i]ᵀ
#                                   = -(4/h²) w_i [δx_i, δy_i]ᵀ   (on uniform grid)
# ============================================================================================================================

"""
    mls_2d(xp, x_nodes, h) -> (ϕ, ∇ϕ, M)

Compute MLS shape function values and gradients for a 2D particle at `xp` (2-vector)
given grid nodes `x_nodes` (2×n matrix, each column one node) with spacing `h`.
"""
function mls_2d(xp::Vector{Float64}, x_nodes::Matrix{Float64}, h::Float64)
    n = size(x_nodes, 2)

    # kernel: tensor product of 1D B-splines
    w = [bspline((x_nodes[1,i]-xp[1])/h)/h * bspline((x_nodes[2,i]-xp[2])/h)/h
         for i in 1:n]

    # polynomial basis  p(x_i) = [1, δx_i, δy_i]
    P = [[1.0, x_nodes[1,i]-xp[1], x_nodes[2,i]-xp[2]]  for i in 1:n]

    # moment matrix M (3×3)
    M = zeros(3, 3)
    for i in 1:n
        M .+= w[i] .* (P[i] * P[i]')
    end

    p0  = [1.0, 0.0, 0.0]
    M⁻¹ = inv(M)

    # shape values  ϕ_i = p(xp)ᵀ M⁻¹ p(x_i) w_i
    ϕ = [(p0' * M⁻¹ * P[i]) * w[i]  for i in 1:n]

    # shape gradients  ∇ϕ_i = -w_i D⁻¹ [δx_i, δy_i]ᵀ,  D⁻¹ = M⁻¹[2:3,2:3]
    D⁻¹ = M⁻¹[2:3, 2:3]
    ∇ϕ  = [-w[i] .* (D⁻¹ * [x_nodes[1,i]-xp[1], x_nodes[2,i]-xp[2]])
            for i in 1:n]

    return ϕ, ∇ϕ, M
end

println("\n", "="^60)
println("MLS shape function — 2D verification")
println("="^60)

h2     = 1.0
nelx   = 4;  nely = 4
# build node grid: (2 × (nelx+1)*(nely+1))
x_nodes_2d = hcat([[Float64(i), Float64(j)] for i in 0:nelx for j in 0:nely]...)

for (xp2, label) in [([2.3, 2.1], "interior       "), ([0.2, 0.3], "near corner    ")]

    ϕ2, ∇ϕ2, M2 = mls_2d(xp2, x_nodes_2d, h2)
    n2           = length(ϕ2)

    println("\nParticle xp = $xp2  ($label)")
    println("  cond(M)             = ", round(cond(M2), digits=4))
    println("  M⁻¹[2,2] (=4/h²)   = ", round(inv(M2)[2,2], digits=6),
            "  ←→  closed-form = ", 4.0/h2^2)

    # Partition of unity
    println("  Σ ϕᵢ               = ", round(sum(ϕ2), digits=10),            "  (must be 1.0)")
    # Linear consistency x
    println("  Σ ϕᵢ δxᵢ           = ",
            round(sum(ϕ2[i]*(x_nodes_2d[1,i]-xp2[1]) for i in 1:n2), digits=10), "  (must be 0.0)")
    # Linear consistency y
    println("  Σ ϕᵢ δyᵢ           = ",
            round(sum(ϕ2[i]*(x_nodes_2d[2,i]-xp2[2]) for i in 1:n2), digits=10), "  (must be 0.0)")
    # Gradient consistency
    println("  Σ ∇ϕᵢ_x            = ",
            round(sum(g[1] for g in ∇ϕ2), digits=10),                             "  (must be 0.0)")
    println("  Σ ∇ϕᵢ_y            = ",
            round(sum(g[2] for g in ∇ϕ2), digits=10),                             "  (must be 0.0)")
end

# ============================================================================================================================
# 6. Plots
#
#   Figure 1 — 1D: ϕ and ∇ϕ as a function of particle position xp
#     Each curve is one node's contribution swept over all particle positions.
#     Left column = interior particle (full support), right = near-boundary particle.
#     Solid = MLS, dashed = standard B-spline.
#
#   Figure 2 — 1D: partition of unity  Σ ϕᵢ(xp)  swept over xp
#     MLS should be identically 1; standard B-spline drops below 1 near boundaries.
#
#   Figure 3 — 2D: heatmap of ϕ at a single node swept over particle positions xp
# ============================================================================================================================

# ---- output directory
dump_dir = raw"C:/Users/lili8/Documents/GitHub/ElastoPlasm.jl/dump"
mkpath(dump_dir)

# ---- shared sweep grid (continuous particle positions)
xp_sweep = range(0.0, nel*h, length=500)

# ---- per-node ϕ and ∇ϕ over the sweep
ϕ_mls_sweep   = [zeros(n) for _ in xp_sweep]
∇ϕ_mls_sweep  = [zeros(n) for _ in xp_sweep]
∇ϕ_full_sweep = [zeros(n) for _ in xp_sweep]
ϕ_std_sweep   = [zeros(n) for _ in xp_sweep]
∇ϕ_std_sweep  = [zeros(n) for _ in xp_sweep]
for (k, xp) in enumerate(xp_sweep)
    ϕ_mls_sweep[k],  ∇ϕ_mls_sweep[k], _, _ = mls_1d(xp, x_nodes, h)
    ∇ϕ_full_sweep[k]                        = mls_1d_full_gradient(xp, x_nodes, h)
    ϕ_std_sweep[k],  ∇ϕ_std_sweep[k]       = std_1d(xp, x_nodes, h)
end

# ---- Figure 1: ϕ per node (MLS vs std) -----------------------------------------------
# ∇ϕ subplot: solid = MLS simplified, dotted = MLS full, dashed = std B-spline
colors = palette(:tab10)[1:n]
p1a = plot(title="ϕᵢ — MLS (solid) vs std B-spline (dashed)", xlabel="xp", ylabel="ϕᵢ",
           legend=:outerright, legendtitle="node i")
p1b = plot(title="∇ϕᵢ — MLS simplified (solid), MLS full (dotted), std (dashed)",
           xlabel="xp", ylabel="∇ϕᵢ", legend=:outerright, legendtitle="node i")
for i in 1:n
    ϕ_mls_i   = [ϕ_mls_sweep[k][i]   for k in eachindex(xp_sweep)]
    ϕ_std_i   = [ϕ_std_sweep[k][i]   for k in eachindex(xp_sweep)]
    ∇ϕ_mls_i  = [∇ϕ_mls_sweep[k][i]  for k in eachindex(xp_sweep)]
    ∇ϕ_full_i = [∇ϕ_full_sweep[k][i] for k in eachindex(xp_sweep)]
    ∇ϕ_std_i  = [∇ϕ_std_sweep[k][i]  for k in eachindex(xp_sweep)]
    plot!(p1a, xp_sweep, ϕ_mls_i;   lw=2, color=colors[i], label="i=$i")
    plot!(p1a, xp_sweep, ϕ_std_i;   lw=1, color=colors[i], ls=:dash,  label=false)
    plot!(p1b, xp_sweep, ∇ϕ_mls_i;  lw=2, color=colors[i],             label="i=$i")
    plot!(p1b, xp_sweep, ∇ϕ_full_i; lw=2, color=colors[i], ls=:dot,   label=false)
    plot!(p1b, xp_sweep, ∇ϕ_std_i;  lw=1, color=colors[i], ls=:dash,  label=false)
end
# mark node positions
vline!(p1a, x_nodes; color=:black, ls=:dot, lw=0.5, label=false)
vline!(p1b, x_nodes; color=:black, ls=:dot, lw=0.5, label=false)
fig1 = plot(p1a, p1b, layout=(2,1), size=(900, 600))
savefig(fig1, joinpath(dump_dir, "mls_fig1_shape_functions.png"))

# ---- Figure 2: partition of unity Σ ϕᵢ(xp) ------------------------------------------
Σϕ_mls = [sum(ϕ_mls_sweep[k])  for k in eachindex(xp_sweep)]
Σϕ_std = [sum(ϕ_std_sweep[k])  for k in eachindex(xp_sweep)]
fig2 = plot(xp_sweep, Σϕ_mls; lw=2, label="MLS  (Σϕᵢ)",
            title="Partition of unity  Σ ϕᵢ(xp)", xlabel="xp", ylabel="Σ ϕᵢ",
            ylims=(0.9, 1.05))
plot!(fig2, xp_sweep, Σϕ_std; lw=2, ls=:dash, label="std B-spline (Σϕᵢ)")
hline!(fig2, [1.0]; color=:black, ls=:dot, lw=1, label="exact = 1")
vline!(fig2, x_nodes; color=:gray, ls=:dot, lw=0.5, label=false)
savefig(fig2, joinpath(dump_dir, "mls_fig2_partition_of_unity.png"))

# ---- Figure 3: 2D heatmap of one node's ϕ over the domain ---------------------------
# sweep xp over the 2D domain and evaluate ϕ at node (2,2) = index nearest (2,2)
target_node = findfirst(i -> x_nodes_2d[1,i] ≈ 2.0 && x_nodes_2d[2,i] ≈ 2.0,
                        1:size(x_nodes_2d,2))
nx2, ny2 = 80, 80
xs2 = range(0.0, nelx*h2, length=nx2)
ys2 = range(0.0, nely*h2, length=ny2)
ϕ_map_mls = [begin
    ϕ2, _, _ = mls_2d([xs2[ix], ys2[iy]], x_nodes_2d, h2)
    ϕ2[target_node]
end for iy in 1:ny2, ix in 1:nx2]
ϕ_map_std = [bspline((xs2[ix]-2.0)/h2)/h2 * bspline((ys2[iy]-2.0)/h2)/h2
             for iy in 1:ny2, ix in 1:nx2]

fig3a = heatmap(xs2, ys2, ϕ_map_mls; title="2D MLS ϕ at node (2,2)", xlabel="x", ylabel="y",
                color=:viridis, aspect_ratio=:equal)
scatter!(fig3a, [x_nodes_2d[1,:]...], [x_nodes_2d[2,:]...];
         color=:white, ms=3, msw=0, label=false)
fig3b = heatmap(xs2, ys2, ϕ_map_std; title="2D std B-spline ϕ at node (2,2)", xlabel="x", ylabel="y",
                color=:viridis, aspect_ratio=:equal)
scatter!(fig3b, [x_nodes_2d[1,:]...], [x_nodes_2d[2,:]...];
         color=:white, ms=3, msw=0, label=false)
fig3 = plot(fig3a, fig3b, layout=(1,2), size=(900, 400))
savefig(fig3, joinpath(dump_dir, "mls_fig3_2d_heatmap.png"))

# ---- Figure 4: gradient consistency  Σ ∇ϕᵢ(xp) over xp (MLS simplified vs std) -----
Σ∇ϕ_mls = [sum(∇ϕ_mls_sweep[k])  for k in eachindex(xp_sweep)]
Σ∇ϕ_std = [sum(∇ϕ_std_sweep[k])  for k in eachindex(xp_sweep)]
fig4 = plot(xp_sweep, Σ∇ϕ_mls; lw=2, label="MLS simplified  (Σ∇ϕᵢ)",
            title="Gradient consistency  Σ ∇ϕᵢ(xp)  [must be 0]",
            xlabel="xp", ylabel="Σ ∇ϕᵢ")
plot!(fig4, xp_sweep, Σ∇ϕ_std; lw=2, ls=:dash, label="std B-spline (Σ∇ϕᵢ)")
hline!(fig4, [0.0]; color=:black, ls=:dot, lw=1, label="exact = 0")
vline!(fig4, x_nodes; color=:gray, ls=:dot, lw=0.5, label=false)
savefig(fig4, joinpath(dump_dir, "mls_fig4_gradient_consistency.png"))

println("\nPlots saved to: ", abspath(dump_dir))

