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
function get_dt(mpts::Point{T1,T2,D,CM},props::MeshProperties{T1,T2,D},time::Time{T1,T2},ΔT::T2) where {T1,T2,D,CM}
    # per-particle elastic wave speed, max over particles (real correctness improvement
    # over the old single global `cmpr.c` scalar if `cmp`/`ρ` ever become heterogeneous)
    cwave = T2(0.0)
    @inbounds for p ∈ 1:mpts.nmp
        cmp   = mpts.s.cmp[p]
        cwave = max(cwave, sqrt((cmp.Kc+T2(4.0/3.0)*cmp.Gc)/mpts.s.ρ[p]))
    end
    # calculte dt
    cmax = props.h./(mpts.vmax.+cwave); mpts.vmax.=T2(0.0)
    dt   = min(T2(0.5)*maximum(cmax),ΔT-time.t[1])
    return dt::T2
end

"""
    get_g(mesh::Mesh{T1,T2,D}; G::T2=9.81) -> Vector{T2}

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
function get_g(props::MeshProperties{T1,T2,D}; G::T2=9.81) where {T1,T2,D}
    if D == 1
        g = [-G]
    elseif D == 2
        g = [T2(0.0),-G]
    elseif D == 3
        g = [T2(0.0),T2(0.0),-G]
    end
    return g::Vector{T2}
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