# `finite_DP` never reset `mpts.s.Δλ[p]` on an elastic step

**Status: fixed.**

A particle that yielded once stayed permanently flagged "yielding" for
`nonlocal.jl`'s gate. Same one-line reset its infinitesimal sibling already
had.
