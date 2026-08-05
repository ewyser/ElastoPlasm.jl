# Pseudo-transient (DYREL/Chebyshev) implicit MPM solver
# Adapted from the FEMTools.jl PT framework (thermal & Stokes DYREL)
#
# Solves the nonlinear momentum residual R = 0 for nodal velocity v^{n+1}:
#   R_i = f_int_i(v^{n+1}) + f_ext_i - m_i * (v_i^{n+1} - v_i^n) / dt
#
# Preconditioner: PC_i = m_i / dt  (diagonal inertia)
# λ_max: from nodal mass spectrum   (Gershgorin bound)
# λ_min: Rayleigh quotient updated every ncheck iterations
# Chebyshev step: α = 2Δτ²/(2 + c·Δτ),  β = (2 − c·Δτ)/(2 + c·Δτ)
#
# NOTE: only LinearElasticity is supported. FiniteElasticity requires
#       additional handling of the bᵢⱼ left-Cauchy-Green tensor.



function pt_solve!(mpts, mesh, cmpr, g, dt, instr;
    iterMax::Int = 500,
    ncheck ::Int = 50,
    CFL         = 0.9,
    c_fact      = 1.0,
    ϵ           = 1e-6)

    T    = eltype(mesh.s.m)
    nno  = mesh.prprt.nno[end]
    ndim = mesh.prprt.dim

    # save t^n state — restored before returning so elasto! can do the permanent update
    mv_n = copy(mesh.s.mv)   # m_i * v_i^n
    σ_n  = copy(mpts.s.σᵢ)   # stress at t^n

    # PT workspace
    ∂v∂τ = zeros(T, ndim, nno)
    R    = zeros(T, ndim, nno)
    R0   = zeros(T, ndim, nno)

    # λ_max of the preconditioned system (PC = M/dt):
    #   precond. Jacobian eigenvalue ≈ 1 + (c_el·dt/h)²  (dimensionless CFL²)
    cfl_mpm  = cmpr.c * dt / minimum(mesh.prprt.h)
    λmax_prc = T(1) + T(cfl_mpm)^2
    Δτ       = T(2) / sqrt(λmax_prc) * T(CFL)  # dimensionless pseudo-step
    α        = Δτ
    β        = T(0)
    nr0  = T(0)

    for it in 1:iterMax
        do_check = (mod(it, ncheck) == 0) || it == 1
        do_check && copyto!(R0, R)

        # momentum residual: R_i = oobf_i - m_i*(v_i^k - v_i^n)/dt
        @inbounds for no in 1:nno
            mi = mesh.s.m[no]
            iszero(mi) && continue
            for d in 1:ndim
                v_n_d = mv_n[d,no] / mi
                r     = mesh.s.oobf[d,no] - mi * (mesh.s.v[d,no] - v_n_d) / dt
                mesh.s.bcs.status[d,no] && (r = T(0); ∂v∂τ[d,no] = T(0))
                R[d,no] = r
            end
        end

        # Chebyshev update: ∂v∂τ ← R/PC + β·∂v∂τ,  v ← v + α·∂v∂τ
        @inbounds for no in 1:nno
            iszero(mesh.s.m[no]) && continue
            pc = mesh.s.m[no] / dt
            for d in 1:ndim
                ∂v∂τ[d,no]     = R[d,no]/pc + β * ∂v∂τ[d,no]
                mesh.s.v[d,no] += α * ∂v∂τ[d,no]
                mesh.s.bcs.status[d,no] && (mesh.s.v[d,no] = T(0))
            end
        end

        # re-evaluate trial stress from updated v^k:
        #   trial_deform: writes ΔFᵢⱼ from current mesh.s.v (no F/J/ρ/Ω update)
        #   reset σ to σ^n, then elast! computes σ^k = σ^n + Δσ(ΔF^k)
        instr.cairn.implicit.trial_deform!(mpts,mesh,dt; ndrange=mpts.nmp);sync(CPU())
        mpts.s.σᵢ .= σ_n
        instr.cairn.implicit.elast!(mpts,cmpr.Del; ndrange=mpts.nmp);sync(CPU())

        # recompute f_int from trial stress (mass is fixed)
        fill!(mesh.s.oobf, T(0))
        instr.cairn.implicit.fint_p2n!(mpts,mesh,g; ndrange=mpts.nmp);sync(CPU())

        # convergence check and λ_min update (Rayleigh quotient)
        if do_check
            nr = norm(R)
            it == 1 && (nr0 = max(nr, eps(T)))
            converged = nr / nr0 < ϵ
            @info "PT" iter=it residual=nr/nr0 α=α β=β
            converged && break

            ΔR    = R .- R0
            denom = sum(abs2, ∂v∂τ) * α^2
            if !iszero(denom)
                numer = T(0)
                @inbounds for no in 1:nno
                    iszero(mesh.s.m[no]) && continue
                    pc = mesh.s.m[no] / dt
                    for d in 1:ndim
                        numer += α * ∂v∂τ[d,no] * ΔR[d,no] / pc
                    end
                end
                λmin = abs(numer) / denom
                c    = T(2) * sqrt(max(λmin, T(0))) * T(c_fact)
                α    = T(2) * Δτ^2 / (T(2) + c * Δτ)
                β    = (T(2) - c * Δτ) / (T(2) + c * Δτ)
            end
        end
    end

    # compute nodal acceleration for FLIP n2p: a = (v^{n+1} - v^n) / dt
    @inbounds for no in 1:nno
        iszero(mesh.s.m[no]) && continue
        for d in 1:ndim
            v_n_d = mv_n[d,no] / mesh.s.m[no]
            mesh.s.a[d,no] = (mesh.s.v[d,no] - v_n_d) / dt
        end
    end

    # restore σ^n so that the subsequent elasto! call does the permanent update
    mpts.s.σᵢ .= σ_n

    return nothing
end

function init_implicit(instr::NamedTuple; implicit::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}())
    implicit[:trial_deform!] = trial_deform(CPU())
    if instr[:perf][:status]
        implicit[:elast!] = elast_fast(CPU())
    else
        implicit[:elast!] = elast(CPU())
    end
    implicit[:fint_p2n!] = fint_p2n(CPU())
    implicit[:n2p!]      = picflip_n2p(CPU())
    return (; implicit...)
end

function mapsto_implicit(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,g::Vector{T2},dt::T2,instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R}
    # reset and accumulate: mass, initial momentum, initial f_int from σ^n
    fill!(mesh.s.m   ,T2(0))
    fill!(mesh.s.mv  ,T2(0))
    fill!(mesh.s.oobf,T2(0))
    fill!(mesh.s.a   ,T2(0))
    fill!(mesh.s.v   ,T2(0))
    instr.cairn.mapsto.map.p2n!(mpts,mesh,g; ndrange=mpts.nmp);sync(CPU())

    # set v = v_n so that the initial PT residual equals the out-of-balance force
    @inbounds for no in 1:mesh.prprt.nno[end]
        iszero(mesh.s.m[no]) && continue
        mi = mesh.s.m[no]
        for d in 1:mesh.prprt.dim
            mesh.s.v[d,no] = mesh.s.mv[d,no] / mi
            mesh.s.bcs.status[d,no] && (mesh.s.v[d,no] = T2(0))
        end
    end

    # PT solve: iterates to find converged v^{n+1}, restores σ^n on exit
    pt_solve!(mpts,mesh,cmpr,g,dt,instr)

    # map converged v^{n+1} back to material points (PIC/FLIP)
    instr.cairn.implicit.n2p!(mpts,mesh,dt,T2(instr.fwrk.C_pf); ndrange=mpts.nmp);sync(CPU())

    return nothing
end
