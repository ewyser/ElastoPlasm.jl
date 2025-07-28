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
function get_dt(mp,mesh,cmpr,instr,time,ΔT)
    # check elastoplastic status
    if time.t[1] > time.te 
        instr[:plast][:status] = true 
    end
    # calculte dt
    cmax = mesh.h./(mp.vmax.+cmpr[:c]); mp.vmax.=0.0 
    dt   = min(0.5*maximum(cmax),ΔT-time.t[1])
    return dt
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
function get_g(dim; G::Float64=9.81)
    if dim == 1 
        g = [-G] 
    elseif dim == 2 
        g = [0.0,-G] 
    elseif dim == 3 
        g = [0.0,0.0,-G] 
    end
    return g
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
function get_spacetime(mp,mesh,cmpr,instr,time,ΔT)
    # calculte dt
    dt = get_dt(mp,mesh,cmpr,instr,time,ΔT)
    # ramp-up gravity
    if time.t[1] <= time.tg 
        g = get_g(mesh.dim; G = 9.81*time.t[1]/time.tg)
    else
        g = get_g(mesh.dim;)
    end
    return g,dt
end