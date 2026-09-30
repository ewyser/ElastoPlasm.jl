# Implicit quasi-static MPM by dynamic relaxation (DR), after Sutcliffe, Coombs, Zhang, Duddu,
# Pretti & Augarde (2026), "A Matrix-Free, Scalable, and Robust Implicit Material Point
# Method", IJNME 127:e70384 — Algorithm 1 (no aggregation). Equation numbers below refer to
# that paper.

function init_implicit(instr::NamedTuple; implicit::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}())
    implicit[:mass_assembly!] = mass_assembly(CPU())
    implicit[:kinematics!]    = dr_kinematics(CPU())
    implicit[:fint!]          = dr_fint(CPU())
    implicit[:update!]        = update(CPU())
    return (; implicit...)
end

"""
    fictitious_mass(mpts, mesh, basis; ζM=2) -> Vector

Diagonal fictitious mass `M̃ᵥᵥ = ζ_M ΔT²/4 Σₚ (U_est/h)² M′ₚ Vₚ` (Eq 35, `ΔT = 1`), summed
over the material points whose support contains node `v`, with the Hencky P-wave modulus
`M′ₚ = exp(−ε_III)/det F · (K + 4G/3)` (Eqs 38, 50). Assembled at the start of the load step,
where `ΔF = I`, so `U_est = 1` (Eq 36). `ε_III` is the smallest principal logarithmic strain;
in 2D plane strain the out-of-plane principal strain is zero and takes part in the minimum.
"""
function fictitious_mass(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D}; ζM::T2=T2(2.0)) where {T1,T2,D}
    h⁻² = inv(minimum(mesh.prprt.h))^2
    M̃   = zeros(T2, mesh.prprt.nno[end])
    for p ∈ 1:mpts.nmp
        λ    = eigen(mpts.s.ϵn[p]).values
        ϵIII = D == 2 ? min(minimum(λ), T2(0.0)) : minimum(λ)
        cmp  = mpts.s.cmp[p]
        M′   = exp(-ϵIII) / mpts.J[p] * (cmp.Kc + T2(4.0/3.0)*cmp.Gc)
        for no ∈ basis.p2n[p]
            iszero(no) && continue
            M̃[no] += h⁻² * M′ * mpts.Ω[p]
        end
    end
    return ζM / T2(4.0) .* M̃
end

"""
    relax!(mpts, mesh, basis, g, instr; tol=1e-9, itmax=1_000_000) -> Int

Solve one quasi-static load step (Algorithm 1) and return the number of DR iterations. The
nodal displacement increment `Δu` is the primary variable and is held in `mesh.s.v`. Each
iteration recomputes stress and out-of-balance force `−R = fᵉˣᵗ − fⁱⁿᵗ` from the current `Δu`,
picks the critical damping factor `α` (Eq 40), and advances the fictitious velocity and `Δu`
with the midpoint rule (Eqs 32-33, `ΔT = 1`; the denominator is `2 + αΔT`, as in Eq 44),
until `‖R‖/‖fᵉˣᵗ‖ < tol` over the free degrees of freedom.
"""
function relax!(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},g::Vector{T2},instr::S; tol::T2=T2(1e-9), itmax::Int=1_000_000) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    # reset grid, lumped mass and external force fᵉˣᵗ = Σ S m g (Eq 22)
    fill!(mesh.s.m, T2(0.0))
    fill!(mesh.s.v, zero(eltype(mesh.s.v)))
    instr.cairn.implicit.mass_assembly!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
    M̃    = fictitious_mass(mpts,mesh,basis)
    free = .!mesh.s.bcs.status .& (M̃ .> T2(0.0))'
    fᵉˣᵗ = g .* mesh.s.m'
    nfᵉˣᵗ = norm(fᵉˣᵗ)
    # midpoint starter v₋½ = −½ΔT M̃⁻¹ fᵉˣᵗ
    v    = [free[d,no] ? -T2(0.5) * fᵉˣᵗ[d,no] / M̃[no] : T2(0.0) for d ∈ 1:D, no ∈ axes(free,2)]
    oobf₀ = similar(mesh.s.oobf)
    α     = T2(0.0)
    for it ∈ 0:itmax
        # stress and out-of-balance force for the current Δu
        fill!(mesh.s.oobf, T2(0.0))
        instr.cairn.implicit.kinematics!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
        instr.cairn.update.elast!(mpts; ndrange=mpts.nmp);sync(CPU())
        instr.cairn.implicit.fint!(mpts,mesh,basis,g; ndrange=mpts.nmp);sync(CPU())
        r = sqrt(sum(abs2, mesh.s.oobf[free]))
        r / nfᵉˣᵗ < tol && return it
        # critical damping α = 2ω₀ from the Rayleigh quotient (Eqs 39-40), fⁱⁿᵗᵢ − fⁱⁿᵗᵢ₋₁ = −(oobfᵢ − oobfᵢ₋₁)
        if it > 0
            num = den = T2(0.0)
            for no ∈ axes(free,2), d ∈ 1:D
                free[d,no] || continue
                num += v[d,no] * (oobf₀[d,no] - mesh.s.oobf[d,no])
                den += v[d,no]^2 * M̃[no]
            end
            α   = num > T2(0.0) && den > T2(0.0) ? min(T2(2.0) * sqrt(num / den), T2(1.9)) : T2(0.0)
        end
        copyto!(oobf₀, mesh.s.oobf)
        # midpoint update of fictitious velocity (Eq 32) and displacement increment (Eq 33)
        c₁, c₂ = (T2(2.0) - α) / (T2(2.0) + α), T2(2.0) / (T2(2.0) + α)
        for no ∈ axes(free,2)
            Δu = MVector(mesh.s.v[no])
            for d ∈ 1:D
                free[d,no] || continue
                v[d,no] = c₁ * v[d,no] + c₂ * mesh.s.oobf[d,no] / M̃[no]
                Δu[d]  += v[d,no]
            end
            mesh.s.v[no] = SVector(Δu)
        end
    end
    @warn "relax!: no convergence after $itmax iterations"
    return itmax
end
