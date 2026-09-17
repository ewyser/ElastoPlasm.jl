# `dynamic_relaxation`'s u-P path throws if exercised

**Status: open.**

`implicit.jl`'s `pt_solve_uP!` calls
`instr.cairn.implicit.fint_p2n!(...)`, a key `init_implicit` never
registers. `elastoquasistatic!` itself doesn't reach this code path, so it
isn't blocked by this; out of scope until someone designs what `fint_p2n!`
should do.
