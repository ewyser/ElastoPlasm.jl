# `vmax` was updated by a racing atomic store in `n2p`

**Status: fixed** (branch `feat-gpu-explicit`).

`picflip_n2p` (`src/home/core/solver/explicit/mapsto/fun/n2p.jl`) kept the largest particle speed per axis, which
`get_dt` reads for the next step's CFL bound, with

```julia
@atom mpts.vmax[d] = max(mpts.vmax[d], abs(v_new[d]))
```

That is an atomic *store* after a non-atomic read, not an atomic maximum: two threads can read the same `vmax[d]`
and the smaller result can overwrite the larger, so with several threads `vmax` could be underestimated and the time
step taken too large. Single-threaded it was exact. On CUDA it could not run at all: Atomix's CUDA extension
implements the atomic store (`Atomix.set!`) as `error("not implemented")`, and CUDA has no floating-point atomic
maximum.

Fix: `n2p` no longer touches `vmax`; `mapsto` reduces the particles' speeds right after it,
`mapreduce(v -> abs.(v), (a, b) -> max.(a, b), mpts.s.v; init=0)`, and copies the result into `mpts.vmax`. `max` is
exact and does not depend on the order, and `get_dt` zeroes `vmax` at the start of each step, so the values are those
of a correct run, also before the first `n2p` (an initial velocity does not enter the first step, as before).
Single-threaded results are bit-identical to before (harness, 10 cases).
