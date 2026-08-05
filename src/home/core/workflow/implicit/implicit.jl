export elastoimplicit!,elastoquasistatic!,elastouP!,elastoquasistaticuP!

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
    nit::Int         = 5000,
    ncheck::Int      = 50,
    CFL              = 0.9,
    c_fact           = 1.0,
    ϵ                = 1e-6,
    quasi_static     = false)

    T    = eltype(mesh.s.m)
    nno  = mesh.prprt.nno[end]
    ndim = mesh.prprt.dim

    # save t^n state — restored before returning so elasto! can do the permanent update
    σ_n  = copy(mpts.s.σᵢ)   # stress at t^n

    # stiffness-based preconditioner scale for quasi-static mode
    h_min = minimum(mesh.prprt.h)
    c_el² = (cmpr.Kc + T(4)/T(3)*cmpr.Gc) / minimum(mpts.s.ρ)

    # PT workspace
    ∂v∂τ = zeros(T, ndim, nno)
    R    = zeros(T, ndim, nno)
    R0   = zeros(T, ndim, nno)

    # λ_max of the preconditioned system (dimensionless, O(1–2))
    # quasi-static: fully stiffness-preconditioned → λmax ≈ 2
    λmax_prc = T(2)
    # quasi-static: factor 1 (not 2) keeps α·λmax < 2 for β=0 (Jury stability)
    Δτ   = T(1)
    α    = Δτ
    β    = T(0)
    nr0  = T(0)

    for it in 1:nit
        do_check = (mod(it, ncheck) == 0) || it == 1
        do_check && copyto!(R0, R)

        # Evaluate trial stress from updated v^k:
        instr.cairn.implicit.trial_deform!(mpts,mesh,dt; ndrange=mpts.nmp);sync(CPU())
        mpts.s.σᵢ .= σ_n
        instr.cairn.implicit.elast!(mpts,cmpr.Del; ndrange=mpts.nmp);sync(CPU())

        # Compute oobf from trial stress
        fill!(mesh.s.oobf, T(0.0))
        instr.cairn.implicit.oobf_assembly!(mpts,mesh,g; ndrange=mpts.nmp);sync(CPU())

        # Residual
        @inbounds for no in 1:nno
            mi = mesh.s.m[no]
            iszero(mi) && continue
            for d in 1:ndim
                r = mesh.s.oobf[d,no]
                mesh.s.bcs.status[d,no] && (r = T(0); ∂v∂τ[d,no] = T(0))
                R[d,no] = r
            end
        end
        # Chebyshev update: ∂v∂τ ← R/PC + β·∂v∂τ,  v ← v + α·∂v∂τ
        @inbounds for no in 1:nno
            iszero(mesh.s.m[no]) && continue
            pc = c_el² * mesh.s.m[no] / h_min^2
            for d in 1:ndim
                ∂v∂τ[d,no]      = R[d,no]/pc + β * ∂v∂τ[d,no]
                mesh.s.v[d,no] += α * ∂v∂τ[d,no]
                if mesh.s.bcs.status[d,no]
                    mesh.s.v[d,no] = T(0)
                end
            end
        end
        # convergence check and λ_min update (Rayleigh quotient)
        if do_check
            nr = norm(R)
            it == 1 && (nr0 = max(nr, eps(T)))
            converged = nr / nr0 < ϵ
            @info "DR pseudo-transient:" iter=it residual=nr/nr0 α=α β=β
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

    # restore σ^n so that the subsequent elasto! call does the permanent update
    mpts.s.σᵢ .= σ_n

    return nothing
end

function init_implicit(instr::NamedTuple; implicit::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}())
    implicit[:assembly!]     = assembly(CPU())
    implicit[:trial_deform!] = trial_deform(CPU())
    implicit[:elast!]        = elast(CPU())
    implicit[:elast_dev!]    = elast_dev(CPU())
    implicit[:oobf_assembly!] = oobf_assembly(CPU())
    implicit[:div_p2n!]      = div_p2n(CPU())
    implicit[:n2p!]          = n2p(CPU())
    return (; implicit...)
end

function relax(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,g::Vector{T2},dt::T2,instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R}
    fill!(mesh.s.m   ,T2(0))
    fill!(mesh.s.mv  ,T2(0))
    fill!(mesh.s.oobf,T2(0))
    fill!(mesh.s.a   ,T2(0))
    fill!(mesh.s.v   ,T2(0))
    instr.cairn.implicit.assembly!(mpts,mesh,g; ndrange=mpts.nmp);sync(CPU())
    pt_solve!(mpts,mesh,cmpr,g,dt,instr; quasi_static=true)
    instr.cairn.implicit.n2p!(mpts,mesh,dt,T2(instr.fwrk.C_pf); ndrange=mpts.nmp);sync(CPU())
    return nothing
end


function elastoquasistatic!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.te)))
    prog = Progress(length(checks);dt=0.5,desc="Solving elastoquasistatic...",barlen=10)
    lstps = collect(range(0, 9.8; length=10))
    for (it,lstp) ∈ enumerate(lstps)

        tic      = time_ns()


        dt = T2(1)
        g  = T2[T2(0), -lstp]

        ignite(mpts,mesh,instr)
        relax(mpts,mesh,cmpr,g,dt,instr)
        elasto(mpts,mesh,cmpr,dt,instr)
        time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)

        savlot(mpts,mesh,time.t[1],instr)
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end

#=
sim = ic_collapse([5, 10], 0.0, 1.0e4, 80.0, 50.0;cli()...)
elastoplasm(sim; workflow=[elastoquasistatic!])
=#





























# ---------------------------------------------------------------------------
# u-P pseudo-transient solver
# Solves coupled momentum + continuity:
#   R_u = f_int(σ') + ∇p·Ω - m*(v-v_n)/dt   (momentum, deviatoric stress)
#   R_p = ∫N(div v + p/K) dΩ                 (continuity, pressure DOF)
# Two independent PT loops: velocity (Δτ_u) and pressure (Δτ_p).
# ---------------------------------------------------------------------------
function pt_solve_uP!(mpts, mesh, cmpr, g, dt, instr;
    iterMax    ::Int  = 500,
    ncheck     ::Int  = 50,
    CFL              = 0.9,
    c_fact           = 1.0,
    ϵ                = 1e-6,
    quasi_static     = false)

    T    = eltype(mesh.s.m)
    nno  = mesh.prprt.nno[end]
    ndim = mesh.prprt.dim
    Kc   = cmpr.Kc
    Gc   = cmpr.Gc
    h_min = minimum(mesh.prprt.h)

    mv_n = copy(mesh.s.mv)
    σ_n  = copy(mpts.s.σᵢ)
    p_n  = copy(mpts.s.p)

    # PT time steps — dimensionless, O(1)
    cfl_u     = quasi_static ? T(0) : cmpr.c * dt / h_min
    λmax_u    = quasi_static ? T(2) : T(1) + T(cfl_u)^2
    # quasi-static: factor 1 (not 2) keeps α·λmax < 2 for β=0 (Jury stability)
    Δτ_u      = (quasi_static ? T(1) : T(2)) / sqrt(λmax_u) * T(CFL)
    # quasi-static: Δτ_p = CFL·K/h²  (stability bound: 2K/h²); dynamic: h²/(K·dt)
    Δτ_p      = quasi_static ? T(CFL) * Kc / h_min^2 : (iszero(dt) ? T(CFL) * h_min / sqrt(Kc) : T(CFL) * h_min^2 / (Kc * dt))

    ∂v∂τ = zeros(T, ndim, nno)
    ∂p∂τ = zeros(T, nno)
    R_u  = zeros(T, ndim, nno)
    R_p  = zeros(T, nno)
    R_u0 = zeros(T, ndim, nno)
    α_u  = Δτ_u;  β_u = T(0)
    nr0  = T(0)

    for it in 1:iterMax
        do_check = (mod(it, ncheck) == 0) || it == 1
        do_check && (copyto!(R_u0, R_u))

        # --- momentum residual ---
        @inbounds for no in 1:nno
            mi = mesh.s.m[no]
            iszero(mi) && continue
            for d in 1:ndim
                r = mesh.s.oobf[d,no]
                quasi_static || (r -= mi * (mesh.s.v[d,no] - mv_n[d,no]/mi) / dt)
                mesh.s.bcs.status[d,no] && (r = T(0); ∂v∂τ[d,no] = T(0))
                R_u[d,no] = r
            end
        end

        # --- continuity residual ---
        @inbounds for no in 1:nno
            iszero(mesh.s.m[no]) && continue
            R_p[no] = -mesh.s.oobp[no]   # oobp = ∫N(divv + p/K)dΩ already assembled
        end

        # --- velocity Chebyshev update ---
        @inbounds for no in 1:nno
            iszero(mesh.s.m[no]) && continue
            pc_u = quasi_static ? Gc * mesh.s.m[no] / h_min^2 : mesh.s.m[no] / dt
            for d in 1:ndim
                ∂v∂τ[d,no]     = R_u[d,no]/pc_u + β_u * ∂v∂τ[d,no]
                mesh.s.v[d,no] += α_u * ∂v∂τ[d,no]
                mesh.s.bcs.status[d,no] && (mesh.s.v[d,no] = T(0))
            end
        end

        # --- pressure update (explicit, no Chebyshev needed for scalar) ---
        @inbounds for no in 1:nno
            iszero(mesh.s.m[no]) && continue
            mesh.s.p[no] += Δτ_p * R_p[no]
        end

        # scatter nodal p back to material points
        T1 = eltype(mesh.prprt.nno)
        @inbounds for p_idx in 1:mpts.nmp
            pp = T(0)
            for nn in 1:mesh.prprt.nn
                no, N, _ = basis(mpts, mesh, T1(p_idx), T1(nn))
                iszero(no) && continue
                pp += N * mesh.s.p[no]
            end
            mpts.s.p[p_idx] = pp
        end

        # --- trial deformation + deviatoric stress reset + elast_dev ---
        instr.cairn.implicit.trial_deform!(mpts,mesh,dt; ndrange=mpts.nmp);sync(CPU())
        mpts.s.σᵢ .= σ_n
        instr.cairn.implicit.elast_dev!(mpts,cmpr.Gc,p_n; ndrange=mpts.nmp);sync(CPU())

        # --- recompute f_int and continuity residual ---
        fill!(mesh.s.oobf, T(0))
        fill!(mesh.s.oobp, T(0))
        instr.cairn.implicit.fint_p2n!(mpts,mesh,g;   ndrange=mpts.nmp);sync(CPU())
        instr.cairn.implicit.div_p2n!( mpts,mesh,Kc,p_n;  ndrange=mpts.nmp);sync(CPU())

        if do_check
            nr_u = norm(R_u);  nr_p = norm(R_p)
            nr   = sqrt(nr_u^2 + nr_p^2)
            it == 1 && (nr0 = max(nr, eps(T)))
            converged = nr / nr0 < ϵ
            @info "u-P PT:" iter=it residual=nr/nr0 R_u=nr_u R_p=nr_p α=α_u
            converged && break

            ΔR_u  = R_u .- R_u0
            denom = sum(abs2, ∂v∂τ) * α_u^2
            if !iszero(denom)
                numer = T(0)
                @inbounds for no in 1:nno
                    iszero(mesh.s.m[no]) && continue
                    # must match the PC used in the velocity update above
                    pc_u = quasi_static ? Gc * mesh.s.m[no] / h_min^2 : mesh.s.m[no] / dt
                    for d in 1:ndim
                        numer += α_u * ∂v∂τ[d,no] * ΔR_u[d,no] / pc_u
                    end
                end
                λmin = abs(numer) / denom
                c    = T(2) * sqrt(max(λmin, T(0))) * T(c_fact)
                α_u  = T(2) * Δτ_u^2 / (T(2) + c * Δτ_u)
                β_u  = (T(2) - c * Δτ_u) / (T(2) + c * Δτ_u)
            end
        end
    end

    @inbounds for no in 1:nno
        iszero(mesh.s.m[no]) && continue
        for d in 1:ndim
            v_n_d = mv_n[d,no] / mesh.s.m[no]
            mesh.s.a[d,no] = (mesh.s.v[d,no] - v_n_d) / dt
        end
    end

    # τᵢ holds the converged u-P total stress so the workflow can restore it after elasto!
    mpts.s.τᵢ .= mpts.s.σᵢ
    # restore σ^n so elasto! can do the permanent kinematic + stress update
    # keep mpts.s.p at converged value so pressure persists between load steps
    mpts.s.σᵢ .= σ_n

    return nothing
end

function mapsto_uP(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,g::Vector{T2},dt::T2,instr::Instruction{T1,T2,D}; quasi_static::Bool=false) where {T1,T2,D,E,R}
    fill!(mesh.s.m   ,T2(0))
    fill!(mesh.s.mv  ,T2(0))
    fill!(mesh.s.oobf,T2(0))
    fill!(mesh.s.oobp,T2(0))
    fill!(mesh.s.a   ,T2(0))
    fill!(mesh.s.v   ,T2(0))
    fill!(mesh.s.p   ,T2(0))
    instr.cairn.mapsto.map.p2n!(mpts,mesh,g; ndrange=mpts.nmp);sync(CPU())
    @inbounds for no in 1:mesh.prprt.nno[end]
        iszero(mesh.s.m[no]) && continue
        mi = mesh.s.m[no]
        for d in 1:mesh.prprt.dim
            mesh.s.v[d,no] = mesh.s.mv[d,no] / mi
            mesh.s.bcs.status[d,no] && (mesh.s.v[d,no] = T2(0))
        end
    end
    # warm-start nodal pressure from carried particle pressures so the first PT
    # scatter doesn't wipe the load-step history
    let vol = zeros(T2, mesh.prprt.nno[end])
        @inbounds for p_idx in T1(1):T1(mpts.nmp)
            for nn in T1(1):T1(mesh.prprt.nn)
                no, N, _ = basis(mpts, mesh, p_idx, nn)
                iszero(no) && continue
                mesh.s.p[no] += N * mpts.Ω[p_idx] * mpts.s.p[p_idx]
                vol[no]       += N * mpts.Ω[p_idx]
            end
        end
        @inbounds for no in 1:mesh.prprt.nno[end]
            iszero(vol[no]) || (mesh.s.p[no] /= vol[no])
        end
    end
    pt_solve_uP!(mpts,mesh,cmpr,g,dt,instr; quasi_static=quasi_static)
    instr.cairn.implicit.n2p!(mpts,mesh,dt,T2(instr.fwrk.C_pf); ndrange=mpts.nmp);sync(CPU())
    return nothing
end

function elastoquasistaticuP!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R}
    it   = T1(0)
    prog = Progress(40; dt=0.5, desc="Solving elasto quasi-static u-P...", barlen=10)
    lstps = collect(range(0, 9.8; length=40))
    for (k, lstp) ∈ enumerate(lstps)
        tic = time_ns()
        dt  = T2(1)
        g   = T2[T2(0), -lstp]

        ignite(mpts, mesh, instr)
        mapsto_uP(mpts, mesh, cmpr, g, dt, instr; quasi_static=true)
        elasto(mpts, mesh, cmpr, dt, instr)
        time.t[1], it, toc = time.t[1]+dt, it+T1(1), (time_ns()-tic)

        savlot(mpts, mesh, time.t[1], instr)
        next!(prog; showvalues=get_vals(mesh, mpts, it))
    end
    finish!(prog)
    return nothing
end

function elastouP!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,D,E,R}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.te)))
    prog = Progress(length(checks);dt=0.5,desc="Solving elasto u-P...",barlen=10)
    for ck ∈ checks
        while time.t[1] < ck
            tic    = time_ns()
            g, dt  = get_spacetime(mpts,mesh,cmpr,time,ck)

            ignite(mpts,mesh,instr)
            mapsto_uP(mpts,mesh,cmpr,g,dt,instr)
            elasto(mpts,mesh,cmpr,dt,instr)
            time.t[1], it, toc = time.t[1]+dt, it+T1(1), (time_ns()-tic)
        end
        savlot(mpts,mesh,time.t[1],instr)
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end