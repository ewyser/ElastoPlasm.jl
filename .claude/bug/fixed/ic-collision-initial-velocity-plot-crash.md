# `ic_collision`'s initial-velocity plot crashes

**Status: resolved by removal** (branch `refactor-init-problem-spec`). `ic_collision` is now a
`collision_case` built by `setup_simulation`, which plots the initial fields (`plot_initial_fields`)
like every other example; the hand-built velocity plot is gone, and `ic_collision` with plotting on
sets up and saves its figure. A velocity plot, if wanted again, should go through
`get_mpts_variable_config()` (`"v"`).

`what_plot_field` errors with `type NamedTuple has no field data`, because
`collision.jl`'s `plot.what` entries are hand-built `NamedTuple`s missing
the `data`/`label`/`unit`/`scale`/`cb` fields `what_plot_field` expects,
instead of going through `get_mpts_variable_config()[name]` like
`slump_problem` does. Pre-existing, not touched by the `ic_collision`
rewrite that fixed its undefined-`kwargser`/stale-pipeline issues.
