export check_case, estimate_dt

# the zone of each candidate point (sub-cell centres, those of `mpts_populate` and of the drawing window,
# `EasyDraw.subcell_centres`; vertical axis fastest; 0 = none), without building a solver or mesh
function case_zones(s::Case)
    l   = s.layout
    axs = EasyDraw.subcell_centres(l.L, l.nel, l.ni)
    pts = Iterators.product(axs[end], axs[1:end-1]...)
    return [zone_at(l.zones, (t[2:end]..., t[1])) for t ∈ pts][:]
end

"""
    estimate_dt(case) -> (dt, c)

The explicit solver's first time step for a case, as `get_dt` computes it at rest: `0.5 h / c`,
with `c` the largest elastic P-wave speed `√((K + 4G/3)/ρ)` over the zones holding points
(`ρ = (1 − n₀)ρ₀`; improved Hencky's tangent bulk modulus `K/n₀` at zero strain) and `h` the smallest
cell size over the axes.
`(Inf, 0.0)` without points.
"""
function estimate_dt(s::Case; zone=case_zones(s))
    used = filter(>(0), unique(zone))
    isempty(used) && return (Inf, 0.0)
    elastic = s.laws.elastic
    function speed(k)
        p(x) = get(get(s.properties, s.layout.zones[k].name, (;)), x, getfield(s.material, x))
        K, G = get_elastic_moduli(p(:E), p(:nu))
        Kt   = elastic == "improved_hencky" ? K/p(:n0) : K
        return sqrt((Kt + 4G/3)/((1 - p(:n0))*p(:rho0)))
    end
    c = maximum(speed, used)
    L, nel = s.layout.L, s.layout.nel
    return 0.5*minimum(L ./ nel)/c, c
end

"""
    check_case(case) -> (; errors, warnings, info)
    check_case(case, solver) -> (; errors, warnings, info)

What to know before running. For a case, its geometry and physics: an error when no zone covers a sub-cell centre (no
material points); warnings for a zone reaching outside the domain or without points (covered by later zones, or
thinner than a sub-cell); the points and elements. With a solver, also how it solves the case: errors
for an option contradicting the workflow (`solution`, or the quasi-static workflows with
`material.elastic = "hypoelastic"`); warnings for non-square cells under dynamic relaxation (slow convergence),
`tg > te` (the gravity ramp is cut at `te`), 32-bit precision with GIMP (fails at the first step, see
`.claude/bug/known/gimpm-32bit-nan-first-step.md`); the load steps, or the estimated first time step and step count
(`estimate_dt`). The drawing window shows a solver's (`draw_case`).
"""
function check_case(s::Case)
    errors, warnings, info = String[], String[], String[]
    zone = case_zones(s)
    L, nel = s.layout.L, s.layout.nel
    count(>(0), zone) == 0 && push!(errors, "No material points: no zone covers a sub-cell centre of the domain.")
    for (k, z) ∈ enumerate(s.layout.zones)
        outside(z.shape, L) && push!(warnings, "Zone \"$(z.name)\" reaches outside the domain: only the part inside gets material points.")
        any(==(k), zone) || push!(warnings, "Zone \"$(z.name)\" has no material points (covered by later zones, or thinner than a sub-cell).")
    end
    push!(info, "$(count(>(0), zone)) material points, $(prod(nel)) elements ($(length(L))D)")
    return (; errors, warnings, info)
end
function check_case(s::Case, solver::AbstractSolver)
    c = check_case(s)
    errors, warnings, info = copy(c.errors), copy(c.warnings), copy(c.info)
    run = solver.run
    qs  = quasistatic(run.workflow)
    L, nel = s.layout.L, s.layout.nel
    (solver.solution == "implicit") != qs && push!(errors, "The option solution = $(solver.solution) contradicts the workflow $(run.workflow).")
    qs && s.laws.elastic == "hypoelastic" &&
        push!(errors, "elastoquasistatic! is finite strain only: material.elastic = hypoelastic is not supported.")
    h = L ./ nel
    qs && maximum(h) > 1.01*minimum(h) && push!(warnings, "Cells are not square: dynamic relaxation converges slowly on elongated cells.")
    !qs && run.tg > run.te && push!(warnings, "tg > te: the gravity ramp is cut at te.")
    startswith(solver.dtype.precision, "32") && solver.basis.which == "gimpm" &&
        push!(warnings, "32-bit precision with GIMP fails at the first step (known bug).")
    zone = case_zones(s)
    nmp  = count(>(0), zone)
    if qs
        push!(info, "$(run.lstps) load steps, each relaxed to equilibrium.")
    elseif nmp > 0
        dt, c = estimate_dt(s; zone)
        T = run.workflow == "elastodynamic!" ? run.te : run.te + run.tep
        push!(info, "Time step ≈ $(round(dt; sigdigits=3)) s (wave speed $(round(c; sigdigits=3)) m/s): ≈ $(ceil(Int, T/dt)) steps for $(T) s.")
    end
    return (; errors, warnings, info)
end
