# `thermal_problem`'s default `bcs.dirichlet` was a hardcoded 2×2 matrix

**Status: fixed.**

`thermal.jl` had `bcs = (; dirichlet = [:fixed :fixed; :fixed :fixed])` — a
literal 2D-only shape. Any 3D `thermal_problem` call threw `BoundsError:
attempt to access 2×2 Matrix{Symbol} at index [3, 1]` in
`get_bc`/`setup_mesh`, before reaching any basis/transfer-scheme code — 3D
thermal was entirely unreachable through the normal entry point (found while
investigating 3D conformity — see
`.claude/bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md`). Fixed
to `fill(:fixed, length(L), 2)`, one `[lower upper]` row per dimension;
verified both 2D and 3D `thermal_problem`→`elastoplasm!` runs succeed with no
`bcs` override needed.
