# `ic_collision`'s initial-velocity plot crashes

**Status: open.** Workaround: run with `plot.status=false` until fixed.

`what_plot_field` errors with `type NamedTuple has no field data`, because
`collision.jl`'s `plot.what` entries are hand-built `NamedTuple`s missing
the `data`/`label`/`unit`/`scale`/`cb` fields `what_plot_field` expects,
instead of going through `get_mpts_variable_config()[name]` like
`slump_problem` does. Pre-existing, not touched by the `ic_collision`
rewrite that fixed its undefined-`kwargser`/stale-pipeline issues.
