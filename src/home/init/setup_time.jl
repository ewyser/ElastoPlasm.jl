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
function setup_time(solver::S; te=0.0,tg=0.0,tep=0.0) where {S<:AbstractSolver{T1,T2,D}} where {T1<:Integer,T2<:Real,D<:AbstractDimension}
    println("T2 = $(T2)")
    return Time{T1,T2}( 
        T2.([0.0,te+tep]), 
        T2(te), 
        if tg > te T2(te) else T2(tg) end, 
        T2(tep),
    )
end