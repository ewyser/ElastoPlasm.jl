"""
    get_dt(mp, mesh, yd, t, time) -> Float64

Computes the adaptive time step for the simulation based on mesh spacing and material point velocities.

# Arguments
- `mp`: Material point data structure, must contain `vmax`.
- `mesh`: Mesh data structure, must contain `h`.
- `yd`: Small offset added to velocity for stability.
- `t`: Current simulation time.
- `time`: Total simulation time.

# Returns
- `Float64`: The computed time step, limited by remaining simulation time.

# Example
```julia
dt = get_dt(mp, mesh, 1e-6, t, time)
```
"""
function get_dt(mp::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,time::NamedTuple,ΔT::T2) where {T1,T2}
    # calculte dt
    cmax = mesh.h./(mp.vmax.+cmpr[:c]); mp.vmax.=T2(0.0) 
    dt   = min(T2(0.5)*maximum(cmax),ΔT-time.t[1])
    return dt::T2
end

"""
    get_g(t::Float64, tg::Float64, ndim::Int64) -> Vector{Float64}

Calculates the gravity vector for the current time, ramping up to full gravity over duration `tg`.

# Arguments
- `t::Float64`: Current simulation time.
- `tg::Float64`: Gravity ramp duration.
- `ndim::Int64`: Number of spatial dimensions.

# Returns
- `Vector{Float64}`: Gravity vector for the current time and dimension.

# Example
```julia
g = get_g(t, tg, 2)
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
    get_spacetime(mp, mesh, cmpr, instr, t, tg, te, time) -> Tuple{Vector{Float64}, Float64}

Updates simulation state, including plasticity status, adaptive time step, and gravity vector.

# Arguments
- `mp`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmpr`: Compression or constitutive model data.
- `instr`: Instruction/configuration dictionary.
- `t`: Current simulation time.
- `tg`: Gravity ramp duration.
- `te`: End time for gravity ramp.
- `time`: Total simulation time.

# Returns
- `(g, dt)`: Tuple containing the gravity vector and computed time step.

# Example
```julia
g, dt = get_spacetime(mp, mesh, cmpr, instr, t, tg, te, time)
```
"""
function get_spacetime(mp::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,time::NamedTuple,ΔT::T2) where {T1,T2}
    # calculte dt
    dt = get_dt(mp,mesh,cmpr,time,ΔT)
    # ramp-up gravity
    if time.t[1] <= time.tg 
        g = get_g(mesh; G = T2(9.81*time.t[1]/time.tg))
    else
        g = get_g(mesh; G = T2(9.81)                  )
    end
    return g,dt
end