"""
    setup_time(ϵ::Type=Float64; te=0.0, tg=0.0, tep=0.0) -> NamedTuple

Set up the time configuration for a simulation, including total, gravity, and elastoplastic time intervals.

# Arguments
- `ϵ::Type=Float64`: (Optional) Numeric type for time values (default: `Float64`).
- `te`: (Optional) End time for the elastodynamic phase.
- `tg`: (Optional) End time for the gravity ramp-up phase.
- `tep`: (Optional) End time for the elastoplastic phase.

# Returns
- `NamedTuple`: Contains time vector `t`, and individual times `te`, `tg`, and `tep`.

# Example
```julia
time = setup_time(Float64; te=10.0, tg=5.0, tep=2.0)
println(time.t)  # [0.0, 12.0]
```
"""
function setup_time(ϵ::Type=Float64; te=0.0,tg=0.0,tep=0.0)
    time  = (; 
        t = ϵ.([0.0,te+tep]), 
        te = ϵ(te), 
        tg = if tg > te ϵ(te) else ϵ(tg) end, 
        tep = ϵ(tep),
    )
    return time
end