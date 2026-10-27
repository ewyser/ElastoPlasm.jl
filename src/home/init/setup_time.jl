"""
    setup_time(solver; te=0.0, tg=0.0, tep=0.0) -> Time

The time stages of a run, in the solver's float type: an elastic stage until `te` with gravity
ramped up over `tg` (cut at `te`), then a plastic stage of `tep`; `time.t = [0, te + tep]` holds the
current and the final time.
"""
function setup_time(solver::S; te=0.0,tg=0.0,tep=0.0) where {T1<:Integer,T2<:Real,D,S<:AbstractSolver{T1,T2,D}}
    return Time{T1,T2}(
        T2.([0.0,te+tep]),
        T2(te),
        if tg > te T2(te) else T2(tg) end,
        T2(tep),
    )
end