# The one-element-wide column with B-splines and the explicit solver goes NaN within 8 steps

**Status: open, not investigated** (may be inherent rather than a bug).

`column_problem([1, 20], 0.0, 1.0e4, 80.0, 50.0)` with the default solver options (`bsmpm`,
`locking = true`, `musl = true`, damping 0.1, DP, Hencky) run by `elastoplasm!(j;
workflows=[elastodynamic!])` throws `InexactError: trunc(Int64, NaN)` in `element_to_nodes_topology`
(`p2e2n`, `ignite`) at about the 8th step: a particle position became NaN. Found while making the 3D
column smoke test of `refactor-init-problem-spec`; not checked against an older commit.

The supported configuration is the one `test_column.jl` uses (GIMP, no locking/MUSL/damping,
`elastoquasistatic!`), which passes 4/4 and gives the same stress profile in 2D and 3D (error 4.58e-3
at `nel = [1, 20]` and `[1, 1, 20]`, `dump/column3d/column_2d_3d.png` at the time). A column one
element wide gives a cubic B-spline basis almost no support across the width, so the failure may be
the basis's limit rather than a defect; deciding needs a run with GIMP + explicit and with B-splines on
a wider column.
