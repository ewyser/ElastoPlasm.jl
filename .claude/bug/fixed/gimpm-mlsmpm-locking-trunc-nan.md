# `basis.which` ∈ {`"gimpm"`, `"mlsmpm"`} combined with `stab.locking=true` crashes with `InexactError: trunc(Int64, NaN)`

**Status: fixed.**

`volumetric.jl`'s `ΔJp` now has the same `iszero(mesh.s.m[no])` guard the
velocity kernels already had.
