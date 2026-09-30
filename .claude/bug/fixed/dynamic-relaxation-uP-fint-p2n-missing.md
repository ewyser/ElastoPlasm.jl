# `dynamic_relaxation`'s u-P path throws if exercised

**Status: fixed** on `feat-implicit-dr-sutcliffe`.

`implicit.jl`'s `pt_solve_uP!` calls
`instr.cairn.implicit.fint_p2n!(...)`, a key `init_implicit` never
registers. `elastoquasistatic!` itself doesn't reach this code path, so it
isn't blocked by this; out of scope until someone designs what `fint_p2n!`
should do.

## Resolution

Removed rather than repaired: `pt_solve_uP!`, `mapsto_uP`, `elastoquasistaticuP!`, `elastouP!`,
`div.jl` and `n2p.jl` were deleted when `pt_solve!` was replaced by the Sutcliffe et al. (2026)
DR solver. `MeshSolidPhase.oobp`/`.p` are now unused (listed under dead code in
`planned-improvements.md`). A mixed u–P formulation, if wanted, starts from scratch.
