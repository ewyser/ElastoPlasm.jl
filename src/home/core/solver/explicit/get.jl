@inline function _Ktan(::Type{<:AbstractElasticLaw}, strain::AbstractStrain, cmp::AbstractConstitutiveModel{T2}, n₀::T2) where {T2}
    return cmp.Kc
end
@inline function _Ktan(::Type{<:ImprovedHenckySolid}, strain::LogarithmicStrain{T2}, cmp::AbstractConstitutiveModel{T2}, n₀::T2) where {T2}
    ϵ = strain.vol
    n = T2(1.0) - (T2(1.0) - n₀) / exp(ϵ)
    return cmp.Kc / (T2(2.0)*n^3) * (n^2*(ϵ^2 + T2(4.0)*ϵ + T2(2.0)) - ϵ*n*(T2(3.0)*ϵ + T2(4.0)) + T2(2.0)*ϵ^2)
end

"""
    get_dt(mpts, mesh, time, ΔT) -> Float64

Compute the adaptive time step for the simulation based on mesh spacing and material point velocities.

# Arguments
- `mpts`: Material point data structure, must contain `vmax` and `s.cmp`/`s.ρ` (per-particle
  constitutive constants + density, used to compute the per-particle elastic wave speed).
- `mesh`: Mesh data structure, must contain `h`.
- `time`: Named tuple with current and phase times.
- `ΔT`: End time for the current phase or time window.

# Returns
- `Float64`: The computed time step, limited by remaining simulation time.

# Example
```julia
dt = get_dt(mpts, mesh, time, ΔT)
```
"""
function get_dt(mpts::Point{T1,T2,D,CM,ST,EL},props::MeshProperties{T1,T2,D},time::Time{T1,T2},ΔT::T2) where {T1,T2,D,CM,ST,EL}
    # per-particle elastic wave speed, max over particles: a reduction, so the arrays may live on any device (over a
    # lazy broadcast: no temporary array, as fast as a loop)
    c = Broadcast.instantiate(Broadcast.broadcasted(mpts.s.cmp, mpts.s.ϵᵢ, mpts.n₀, mpts.s.ρ) do cmp, ϵ, n₀, ρ
        sqrt((_Ktan(EL,ϵ,cmp,n₀)+T2(4.0/3.0)*cmp.Gc)/ρ)
    end)
    cwave = reduce(max, c; init=T2(0.0))
    # dt: the CFL bound of the most restrictive axis, Courant number 0.5 (Nakamura et al. 2023, Eq 35); vmax, the
    # largest speed per axis (mapsto), read once and zeroed for the next step
    vmax = SVector{D,T2}(Array(mpts.vmax)); fill!(mpts.vmax, T2(0.0))
    cmax = props.h./(vmax.+cwave)
    dt   = min(T2(0.5)*minimum(cmax),ΔT-time.t[1])
    return dt::T2
end

"""
    get_g(props::MeshProperties{T1,T2,D}; G::T2=9.81) -> SVector{D,T2}

Calculate the gravity vector for the mesh, with magnitude `G` in the negative last direction.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object containing dimension information.
- `G::T2=9.81`: (Optional) Gravity magnitude (default: 9.81).

# Returns
- `Vector{T2}`: Gravity vector for the mesh dimension.

# Example
```julia
g = get_g(mesh)
```
"""
function get_g(props::MeshProperties{T1,T2,D}; G::T2=T2(9.81)) where {T1,T2,D}
    return SVector{D,T2}(ntuple(d -> d == D ? -G : T2(0.0), Val(D)))   # a plain value: a kernel argument on any device
end

"""
    get_spacetime(mpts, mesh, time, ΔT) -> Tuple{Vector, Float64}

Update simulation state, including plasticity status, adaptive time step, and gravity vector, for the current time step.

# Arguments
- `mpts`: Material point data structure.
- `mesh`: Mesh data structure.
- `time`: Named tuple with current and phase times.
- `ΔT`: End time for the current phase or time window.

# Returns
- `(g, dt)`: Tuple containing the gravity vector and computed time step for the current time.

# Example
```julia
g, dt = get_spacetime(mpts, mesh, time, ΔT)
```

# Notes
- Computes the adaptive time step using `get_dt`.
- Ramps up gravity linearly until the end of the gravity phase, then applies full gravity.
"""
function get_spacetime(mpts::Point{T1,T2,D,CM},mesh::Mesh{T1,T2,D},time::Time{T1,T2},ΔT::T2) where {T1,T2,D,CM}
    # calculte dt
    dt = get_dt(mpts,mesh.prprt,time,ΔT)
    # ramp-up gravity
    if time.t[1] <= time.tg 
        g = get_g(mesh.prprt; G = T2(9.81*time.t[1]/time.tg))
    else
        g = get_g(mesh.prprt; G = T2(9.81)                  )
    end
    return g,dt
end

"""
    checkpoints(time, solver, stop) -> Vector{Tuple{T2, Bool, Bool}}

The stops of an explicit workflow from the current time to `stop`, as `(t, plot, record)`: the
`plot.freq` grid (a live plot) and, with `output.status`, the `output.every` grid (an output frame),
each ending at `stop` itself when the grid misses it, so the run ends at `stop` rather than at the last
grid point before it (with `stop < plot.freq` the grid holds only the start: no step at all). Without
output the stops are the plot grid alone.
"""
function checkpoints(time::Time{T1,T2}, solver, stop) where {T1,T2}
    grid(step) = (ts = collect(time.t[1]:step:stop); isempty(ts) || last(ts) ≈ stop || push!(ts, stop); T2.(ts))
    plots = grid(solver.plot.freq)
    solver.output.status || return [(t, true, false) for t ∈ plots]
    recs  = grid(solver.output.every)
    return [(t, t ∈ plots, t ∈ recs) for t ∈ sort(unique([plots; recs]))]
end
