export Time

"""
    Time{T1,T2}

The time stages of a run: an elastic stage until `te`, gravity ramped up over its first `tg`, then a plastic stage
of `tep`.

# Fields
- `t::Vector{T2}`: the current time, then the final one, `te + tep` [s].
- `te::T2`: the elastic stage's end [s].
- `tg::T2`: the gravity ramp's duration, at most `te` [s].
- `tep::T2`: the plastic stage's duration [s].

`T1` and `T2` are the index and float types (`T1` names the problem's, no field holds it).

# Constructors
    Time{T1,T2}(t, te, tg, tep)

`setup_time` builds it from a solver's run settings.
"""
struct Time{T1,T2}
    t  ::Vector{T2}
    te ::T2
    tg ::T2
    tep::T2
end
