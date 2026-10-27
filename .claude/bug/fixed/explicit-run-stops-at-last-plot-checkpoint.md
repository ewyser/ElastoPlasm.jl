# Explicit workflows stop at the last `plot.freq` multiple below the end time, not at it

**Status: fixed** (branch `fix-explicit-end-time-and-cfl-step`).

`elastodynamic!`, `elastoplastic!` and `thermodynamic!` (`src/home/core/solver/explicit/worflow.jl`)
step towards the checkpoints `time.t[1]:plot.freq:te` (`elastodynamic!`) or `…:time.t[2]` (the other
two), and nothing else. A range stops at its last element ≤ the end, so the end time itself is a stop
only when it is a multiple of `plot.freq` past the start. Otherwise the run ends short of it; with
`te < plot.freq` (default `freq = 1.0`) the only checkpoint is the start, so the workflow takes **zero
steps** and returns as if it had run. Since `a6bdd9a1` (2026-08-10).

Found when `test_case.jl`'s `execute`/`execute!` tests (a 4×2 block, `te = 0.2`) showed the state
unchanged after a run, with an initial velocity too: the old check `time.t[end] == 0.2` passed only
because `time.t[end]` is `te` from the start. The test now sets `plot.freq = 0.1`. The default
examples (`te = 10`, `tep = 5`, `freq = 1`) are not affected.

Fix: append the end time to the stops (`unique(sort([collect(t:freq:te); te]))`, the plot taken there
too). `feat-hdf5-output`'s `checkpoints` (`explicit/get.jl`) builds the same ranges (`plot.freq` and
`output.every` grids up to `stop`) and needs it there, which is the natural place once that branch is
merged. Changes the end state of any run whose end time is not on the grid, so check that against a
run with a `freq` that divides it.

**Fix:** `checkpoints(time, solver, stop)` (`explicit/get.jl`) gives the three workflows their stops: the
`plot.freq` grid, plus `stop` when the grid misses it (`last ≈ stop`, so an on-grid end adds no stop; the
plot is then also taken at the end). Same name and signature as `feat-hdf5-output`'s, which needs the
same append at the merge. Checked: the default slump (te = 10, tep = 5, freq = 1) bit-identical to the
pre-fix baseline, single-threaded; a slump with te = 1.25, tep = 0.5 now ends its stages at 1.25 s and
1.75 s (before: 1.0 s, then no plastic step); `test_case` asserts that te = 0.2 < freq = 1 runs to 0.2.
