# `gimpm`+`tpic`+`locking=true`+`musl=true` fails only when `nonloc.status=true`

**Status: open, narrow, not root-caused.**

Found as part of the `test_workflow.jl` sweep restructuring documented in
`.claude/bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md`:
`gimpm, finite, tpic, locking=true, musl=true` passes with
`nonloc.status=false` but fails with `nonloc.status=true` (`BoundsError:
attempt to access 400-element Vector{...} at index [16346]` in
`element_to_nodes_topology`, `tplgy.jl:15` — a particle drifting far outside
the mesh, same failure mode as the broader grid-crossing instability, just
newly triggered on an otherwise-stable configuration). Flagged here as a
genuine, narrow, nonlocal-specific instability for `gimpm`+`tpic` worth
investigating alongside the broader grid-crossing theory, since nonlocal
regularization's plastic-strain averaging could plausibly be feeding back
into the same kind of drift.
