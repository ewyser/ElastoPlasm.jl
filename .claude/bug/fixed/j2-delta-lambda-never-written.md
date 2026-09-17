# `mpts.s.Δλ` never written under `plast.constitutive="J2"`, silently disabling non-local regularization

**Status: fixed.**

`J2.jl`'s return-map kernels computed `Δλ` as a purely loop-local quantity
and never assigned `mpts.s.Δλ[p]`; since `nonlocal.jl` gates all three of its
branches on `mpts.s.Δλ[p] != 0`, non-local regularization never activated
under J2. Fix: accumulate `Δλ` across the CPA iterations and store the sum.
Measured before→after on a 591-particle run: `Δλ != 0` on 0→424 particles;
DP runs were never affected.
