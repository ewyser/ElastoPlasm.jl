
"""
    build_solid_phase(T1, T2, D, solver, pts::PointData, v0) -> PointSolidPhase

Build the per-particle solid-phase state (`PointSolidPhase`) of the points `pts` with initial
volumes `v0`: the per-particle constitutive-model bundle (`cmp`, via `setup_cmp`) and all mechanical
state fields, zero-initialized except `v` (`pts.vp`), `ρ₀`, `ρ = (1 − n₀)ρ₀` and `m = ρ·Ω₀` (the
constant solid mass).

`ST` (the typed strain storage of `ϵᵢ`/`ϵn`, see `tensor.jl`) and `EL`
(`Point`'s elastic law, see `AbstractElasticLaw` in `lagrangian.jl`) are
both picked here, together, purely from `solver.material.elastic` — there is no
separate strain-formulation config key. `"hypoelastic"` → `InfinitesimalStrain`+
`HypoelasticSolid`; `"hencky"`/`"improved_hencky"` → `LogarithmicStrain`+
(`HenckySolid`/`ImprovedHenckySolid`). `CM` (the constitutive-model type
of `cmp`) is picked by `solver.material.plastic`: `DruckerPrager` for `"DP"`,
`VonMises` for `"VM"` — see `setup_cmp`.
"""
function build_solid_phase(T1,T2,D,solver,pts::PointData,v0)
    L  = D*D
    nmp, ρ0, vp = length(pts.zone), pts.rho0, pts.vp
    ρ  = T2.((1.0 .- pts.n0) .* ρ0)
    TM = SMatrix{D,D,T2,L}
    TV = SVector{D,T2}
    ST,EL = if solver.material.elastic == "hypoelastic"
        InfinitesimalStrain{T2}, HypoelasticSolid
    elseif solver.material.elastic == "hencky"
        LogarithmicStrain{T2}, HenckySolid
    elseif solver.material.elastic == "improved_hencky"
        LogarithmicStrain{T2}, ImprovedHenckySolid
    else
        throw(error("InvalidElasticLaw: $(solver.material.elastic) (expected \"hypoelastic\", \"hencky\" or \"improved_hencky\")"))
    end
    cmp = setup_cmp(nmp, T2.(pts.c0), T2.(pts.cr), T2.(pts.phi); E=T2.(pts.E), ν=T2.(pts.nu), Hp=T2.(pts.Hp), D=Int(D),
                    constitutive=solver.material.plastic)
    CM  = eltype(cmp)

    s = PointSolidPhase{T1,T2,D,CM,ST,EL,L}(
        [zero(TV)  for _ in 1:nmp]                         , # u
        [TV(T2.(vp[:,p])) for p in 1:nmp]                  , # v
        # mechanical properties
        T2.(ρ0)                                            , # ρ₀
        ρ                                                  , # ρ
        ρ .* T2.(vec(copy(v0)))                            , # m = ρ·Ω₀
        T2.(zeros(nmp))                                    , # Δλ
        [zero(SVector{2,T2}) for _ in 1:nmp]               , # ϵpII
        T2.(zeros(nmp))                                    , # ϵpV
        # typed stress tensors (tensor.jl)
        [zero(CauchyStress{T2}) for _ in 1:nmp]           , # σᵢ
        [zero(CauchyStress{T2}) for _ in 1:nmp]           , # σn
        [zero(KirchhoffStress{T2}) for _ in 1:nmp]        , # τᵢ
        # tensor in matrix notation
        [zero(TM) for _ in 1:nmp]                          , # ∇vᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ∇uᵢⱼ
        [zero(TM) for _ in 1:nmp]                          , # ΔFᵢⱼ
        [TM(I)    for _ in 1:nmp]                          , # Fᵢⱼ
        [TM(I)    for _ in 1:nmp]                          , # Fn
        [zero(ST) for _ in 1:nmp]                          , # ϵᵢ
        [zero(ST) for _ in 1:nmp]                          , # ϵn
        [zero(TM) for _ in 1:nmp]                          , # ωᵢⱼ
        # per-particle constitutive-model constants
        cmp                                                  , # cmp::Vector{CM}
    )
    return s
end

"""
    build_thermal_phase(T1, T2, D, pts::PointData; thermal::Bool=false) -> Union{Nothing,PointThermalPhase}

Build the per-particle thermal-phase state (`PointThermalPhase`) from `pts.cp`/`pts.k`/`pts.T`
when `thermal=true`, otherwise `nothing`.
"""
function build_thermal_phase(T1,T2,D,pts::PointData; thermal::Bool=false)
    return thermal ? PointThermalPhase{T1,T2,D}(
        T2.(pts.cp)                                       , # c::Vector{T2} specific heat capacity vector
        T2.(pts.k)                                        , # k::Vector{T2} thermal conductivity vector
        T2.(zeros(D,length(pts.zone)))                    , # q::Matrix{T2} heat flux array
        T2.(pts.T)                                        , # T::Vector{T2} temperature vector
    ) : nothing
end

"""
    setup_mpts(mesh, solver, pts::PointData; thermal=false) -> Point

Construct the material point system of a simulation from its points (`case_points`):
positions, domains (`ℓ₀`, from the element size and `pts.ni`), volumes `Ω₀`, porosity `n₀`, the
solid phase (`build_solid_phase`) and, with `thermal=true`, the thermal phase. Connectivity is built
separately by `setup_basis`, once both `Mesh` and `Point` exist.
"""
function setup_mpts(mesh::Mesh{T1,T2,D},solver::S,pts::PointData; thermal::Bool=false) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    props   = mesh.prprt
    ni,xp   = pts.ni,pts.xp
    nmp     = size(xp,2)
    l0      = ones(size(xp)).*0.5.*(props.h./ni)
    v0      = prod(2 .* l0; dims=1)
    s       = build_solid_phase(T1,T2,D,solver,pts,v0)
    t       = build_thermal_phase(T1,T2,D,pts; thermal=thermal)
    n0      = pts.n0
    mpts = Point(                              # type parameters inferred from `s`
        # general information
        T1(D)                              , # ndim
        T1(nmp)                              , # nmp
        T2.(zeros(D))                      , # vmax
        # material point properties
        [SVector{D,T2}(xp[:,i]) for i in 1:nmp]           , # x
        [SVector{D,T2}(l0[:,i]) for i in 1:nmp]     , # ℓ₀
        [SVector{D,T2}(l0[:,i]) for i in 1:nmp]     , # ℓ
        T2.(copy(n0))                        , # n₀
        T2.(copy(n0))                        , # n
        T2.(vec(copy(v0)))                   , # Ω₀
        T2.(vec(copy(v0)))                   , # Ω
        T2.(ones(nmp))                       , # ΔJ
        T2.(ones(nmp))                       , # J
        # solid phase
        s                                    , #
        # fluid phase
        nothing                              , #
        # thermal phase
        t                                    , #
    )
    return mpts
end
