# `test_collapse.jl` fully non-functional

**Status: fixed, and since renamed.**

`collapse_problem` (renamed from `ic_collapse`) was rewritten to the current
`get_solver → setup_geometry → setup_problem → setup_basis → export_problem`
pipeline; `test_collapse.jl` updated to read
`file["ic/problem"]`/`file["cfg/solver"]` instead of the old
`ic["mesh"]`/`ic["mpts"]`/`ic["cmpr"]`/`cfg["instr"]` layout. All 4 `nel`
convergence cases pass, including `@test errors[k+1] < errors[k]` — a real
physics-correctness signal. Also the first real exercise of
`dynamic_relaxation`'s `cmp`-rerouted reads (see
`.claude/docs/planned-improvements.md`'s "Typed constitutive-model
abstraction" section), now confirmed working.

**Since renamed**: this test is a 1-D elastic self-weight column convergence
check (no plasticity ever exercised), not a granular collapse — the
misleading name was freed up for an actual granular-collapse example,
`collapse_problem` (see the `collapse_problem`/`column_problem` bullets under
`.claude/docs/operations.md`). `collapse_problem`→`column_problem`,
`collapse.jl`→`column.jl`, `get_collapse`→`get_column`,
`test_collapse.jl`→`test_column.jl`; the 4/4 convergence result above is
unchanged post-rename.
