"""
    get_dt(mpts, mesh, cmpr, time, ΔT) -> Float64

Compute the adaptive time step for the simulation based on mesh spacing and material point velocities.

# Arguments
- `mpts`: Material point data structure, must contain `vmax`.
- `mesh`: Mesh data structure, must contain `h`.
- `cmpr`: Constitutive model parameters (must include wave speed `c`).
- `time`: Named tuple with current and phase times.
- `ΔT`: End time for the current phase or time window.

# Returns
- `Float64`: The computed time step, limited by remaining simulation time.

# Example
```julia
dt = get_dt(mpts, mesh, cmpr, time, ΔT)
```
"""
function get_dt(mpts::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,time::NamedTuple,ΔT::T2) where {T1,T2}
    # calculte dt
    cmax = mesh.h./(mpts.vmax.+cmpr[:c]); mpts.vmax.=T2(0.0) 
    dt   = min(T2(0.5)*maximum(cmax),ΔT-time.t[1])
    return dt::T2
end

"""
    get_g(mesh::Mesh{T1,T2}; G::T2=9.81) -> Vector{T2}

Calculate the gravity vector for the mesh, with magnitude `G` in the negative last direction.

# Arguments
- `mesh::Mesh{T1,T2}`: Mesh object containing dimension information.
- `G::T2=9.81`: (Optional) Gravity magnitude (default: 9.81).

# Returns
- `Vector{T2}`: Gravity vector for the mesh dimension.

# Example
```julia
g = get_g(mesh)
```
"""
function get_g(mesh::Mesh{T1,T2}; G::T2=9.81) where {T1,T2}
    if mesh.dim == T1(1) 
        g = [-G] 
    elseif mesh.dim == T1(2) 
        g = [T2(0.0),-G] 
    elseif mesh.dim == T1(3) 
        g = [T2(0.0),T2(0.0),-G] 
    end
    return g::Vector{T2}
end

"""
    get_spacetime(mpts, mesh, cmpr, time, ΔT) -> Tuple{Vector, Float64}

Update simulation state, including plasticity status, adaptive time step, and gravity vector, for the current time step.

# Arguments
- `mpts`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmpr`: Constitutive model parameters.
- `time`: Named tuple with current and phase times.
- `ΔT`: End time for the current phase or time window.

# Returns
- `(g, dt)`: Tuple containing the gravity vector and computed time step for the current time.

# Example
```julia
g, dt = get_spacetime(mpts, mesh, cmpr, time, ΔT)
```

# Notes
- Computes the adaptive time step using `get_dt`.
- Ramps up gravity linearly until the end of the gravity phase, then applies full gravity.
"""
function get_spacetime(mpts::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,time::NamedTuple,ΔT::T2) where {T1,T2}
    # calculte dt
    dt = get_dt(mpts,mesh,cmpr,time,ΔT)
    # ramp-up gravity
    if time.t[1] <= time.tg 
        g = get_g(mesh; G = T2(9.81*time.t[1]/time.tg))
    else
        g = get_g(mesh; G = T2(9.81)                  )
    end
    return g,dt
end