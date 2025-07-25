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
function get_dt(mp, mesh, yd, t, time)
    h    = mesh.h
    vmax = mp.vmax .+ yd
    if length(h) == length(vmax)
        cmax     = h ./ vmax
        mp.vmax .= 0.0
        return min(0.5 * maximum(cmax), time - t)
    else
        return nothing
    end
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
function get_g(t::Float64,tg::Float64,ndim::Int64)
    g = 0.0
    if t<=tg 
        g = 9.81*t/tg 
    else
        g = 9.81
    end
    return if ndim == 1 g = [-g] elseif ndim == 2 g = [0.0 -g] elseif ndim == 3 g = [0.0 0.0 -g] end
end

"""
    get_spacetime(mp, mesh, cmp, instr, t, tg, te, time) -> Tuple{Vector{Float64}, Float64}

Updates simulation state, including plasticity status, adaptive time step, and gravity vector.

# Arguments
- `mp`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmp`: Compression or constitutive model data.
- `instr`: Instruction/configuration dictionary.
- `t`: Current simulation time.
- `tg`: Gravity ramp duration.
- `te`: End time for gravity ramp.
- `time`: Total simulation time.

# Returns
- `(g, dt)`: Tuple containing the gravity vector and computed time step.

# Example
```julia
g, dt = get_spacetime(mp, mesh, cmp, instr, t, tg, te, time)
```
"""
function get_spacetime(mp,mesh,cmp,instr,t,tg,te,time) # t = tw 
    # check elastoplastic status
    if t > te 
        instr[:plast][:status] = true 
    end
    # calculte dt
    cmax = mesh.h./(mp.vmax.+cmp[:c]); mp.vmax.=0.0 
    dt   = min(0.5*maximum(cmax),time-t)
    # ramp-up gravity
    if t<=tg 
        g = 9.81*t/tg 
    else
        g = 9.81
    end
    if mesh.dim == 1 
        g = [-g] 
    elseif mesh.dim == 2 
        g = [0.0 -g] 
    elseif mesh.dim == 3 
        g = [0.0 0.0 -g] 
    end
    return g,dt
end