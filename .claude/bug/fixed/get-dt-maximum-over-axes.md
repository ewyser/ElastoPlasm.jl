# `get_dt` takes the largest cell size over the axes, not the smallest

**Status: fixed** (branch `fix-explicit-end-time-and-cfl-step`).

`get_dt` (`src/home/core/solver/explicit/get.jl`) computes `cmax = props.h ./ (mpts.vmax .+ cwave)` per
axis and took `dt = 0.5*maximum(cmax)`: the step of the *least* restrictive axis. The CFL bound is the
most restrictive one; Nakamura et al. (2023, *CMAME* 404:115720, Eq 35, `refs/ElastoPlasm/`) give
`Δt_cr = h / max_p(c_p + ‖v_p‖)` and `Δt = Cr Δt_cr`, `Cr ≤ 1`. Found while writing `estimate_dt`, which
mirrored it. Fix: `minimum` in both (`estimate_dt` in `init/case/checks.jl`, still tested equal to
`get_dt` at rest, now also on 0.5 × 1 m cells, where it must be `0.5·0.5/c`).

**Not only elongated cells.** The bug file first said square cells were unaffected. That holds at rest
only: `mpts.vmax` is the peak velocity *per axis* (`n2p.jl`), so `h[d]/(vmax[d] + c)` differs between
axes as soon as points move, and `maximum` took the axis with the slower points.

**Effect** (single thread; the old reduction swapped back in the same session reproduced the pre-fix
baseline bit for bit; figure `dump/fix_dt/dt_minimum_vs_maximum.png`, from a scratch comparison not
kept in the repo, gone once `dump/` is emptied):
- default slump, square 1.604 m cells, 10 s + 5 s: 473 → 476 steps; the same failure (max displacement
  9.557 → 9.532 m, max `ϵpII` 1.877 → 1.902), point positions differing by at most 0.10 m (mean 1.7 cm),
  at the scarp and along the shear band;
- slump on 0.401 × 1.604 m cells, 3 s elastic: 92 → 360 steps (×4, Courant number 2 → 0.5 along x).
  The old step did **not** visibly go unstable in those 3 s (both runs finite, max |v| 0.72 vs 0.70 m/s),
  so the fix rests on the CFL bound, not on an observed blow-up.

`test_workflow`: 145 → 144/192; the one new failure is a gimpm case whose three siblings already failed
(`.claude/bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md`). `test_case` 74/74,
`test_column` 4/4, `test_draw` 105/105.

Still different from Eq 35, left as is: the velocity is the largest component per axis, not the speed
`‖v_p‖`, and the largest wave speed and the largest velocity are taken over the points separately rather
than `max_p(c_p + ‖v_p‖)` (see `planned-improvements.md`).
