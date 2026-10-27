export check_case, estimate_dt

# the zone of each candidate point (sub-cell centres, as `mpts_populate`, vertical axis fastest; 0 =
# none), without building a solver or mesh; an extruded problem is counted in its x–z plane
function case_zones(s::Case)
    L, nel = s.L, s.nel
    axs = [collect(0.5*L[d]/(nel[d]*s.ni):L[d]/(nel[d]*s.ni):L[d] - 0.5*L[d]/(nel[d]*s.ni)) for d ∈ eachindex(L)]
    pts = Iterators.product(axs[end], axs[1:end-1]...)
    return [something(findlast(z -> inside(z.shape, (t[2:end]..., t[1])), s.zones), 0) for t ∈ pts][:]
end
# whether a shape reaches outside the domain [0, L] (the x–z plane for two-coordinate shapes)
function outside(shape::AbstractShape, L)
    lims(c) = length(c) == length(L) ? L : [L[1], L[end]]
    shape isa Polygon && return any(v -> !(0 ≤ v[1] ≤ L[1] && 0 ≤ v[2] ≤ L[end]), shape.vertices)
    shape isa Block     && return any(shape.lo .< 0) || any(shape.hi .> lims(shape.lo))
    shape isa Ball    && return any(shape.centre .- shape.radius .< 0) || any(shape.centre .+ shape.radius .> lims(shape.centre))
    return false
end
# a solver option of the case, or its default
option_of(s::Case, sec, key) = (o = get(s.solver, sec, (;)); o isa NamedTuple && haskey(o, key) ? getfield(o, key) : getfield(getfield(get_default(), sec), key))

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
    elastic = option_of(s, :material, :elastic)
    function speed(k)
        p(x) = get(s.zones[k].props, x, getfield(s.material, x))
        K, G = get_elastic_moduli(p(:E), p(:nu))
        Kt   = elastic == "improved_hencky" ? K/p(:n0) : K
        return sqrt((Kt + 4G/3)/((1 - p(:n0))*p(:rho0)))
    end
    c = maximum(speed, used)
    L, nel = domain(s)
    return 0.5*minimum(L ./ nel)/c, c
end

"""
    check_case(case) -> (; errors, warnings, info)

What to know about a case before running it. Errors: no material points, an option contradicting
the workflow (`solution`, or the quasi-static workflows with `material.elastic = "hypoelastic"`).
Warnings: a zone reaching outside the domain, a zone without points (covered by later zones, or
thinner than a sub-cell), non-square cells under dynamic relaxation (slow convergence), `tg > te`
(the gravity ramp is cut at `te`), 32-bit precision with GIMP (fails at the first step, see
`.claude/bug/known/gimpm-32bit-nan-first-step.md`). Info: points, elements, and the estimated first
time step and step count (`estimate_dt`).
"""
function check_case(s::Case)
    errors, warnings, info = String[], String[], String[]
    zone = case_zones(s)
    run  = s.run
    qs   = quasistatic(run.workflow)
    sol  = get(s.solver, :solution, nothing)
    L, nel = domain(s)
    count(>(0), zone) == 0 && push!(errors, "No material points: no zone covers a sub-cell centre of the domain.")
    !isnothing(sol) && (sol == "implicit") != qs && push!(errors, "The option solution = $sol contradicts the workflow $(run.workflow).")
    qs && option_of(s, :material, :elastic) == "hypoelastic" &&
        push!(errors, "elastoquasistatic! is finite strain only: material.elastic = hypoelastic is not supported.")
    for (k, z) ∈ enumerate(s.zones)
        outside(z.shape, s.L) && push!(warnings, "Zone \"$(z.name)\" reaches outside the domain: only the part inside gets material points.")
        any(==(k), zone) || push!(warnings, "Zone \"$(z.name)\" has no material points (covered by later zones, or thinner than a sub-cell).")
    end
    h = L ./ nel
    qs && maximum(h) > 1.01*minimum(h) && push!(warnings, "Cells are not square: dynamic relaxation converges slowly on elongated cells.")
    !qs && run.tg > run.te && push!(warnings, "tg > te: the gravity ramp is cut at te.")
    startswith(option_of(s, :dtype, :precision), "32") && option_of(s, :basis, :which) == "gimpm" &&
        push!(warnings, "32-bit precision with GIMP fails at the first step (known bug).")
    nmp = count(>(0), zone)*(isnothing(s.width) ? 1 : s.nely*s.ni)
    push!(info, "$nmp material points, $(prod(nel)) elements ($(length(L))D" * (isnothing(s.width) ? ")" : ", extruded along y)"))
    if qs
        push!(info, "$(run.lstps) load steps, each relaxed to equilibrium.")
    elseif nmp > 0
        dt, c = estimate_dt(s; zone)
        T = run.workflow == "elastodynamic!" ? run.te : run.te + run.tep
        push!(info, "Time step ≈ $(round(dt; sigdigits=3)) s (wave speed $(round(c; sigdigits=3)) m/s): ≈ $(ceil(Int, T/dt)) steps for $(T) s.")
    end
    return (; errors, warnings, info)
end
