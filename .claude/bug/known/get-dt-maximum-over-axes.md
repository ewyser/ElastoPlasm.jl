# `get_dt` takes the largest cell size over the axes, not the smallest

**Status: open, not investigated.**

`get_dt` (`src/home/core/solver/explicit/get.jl`) computes `cmax = props.h ./ (mpts.vmax .+ cwave)` per
axis and then `dt = 0.5*maximum(cmax)`. A CFL bound would take the smallest cell size (`minimum`). With square
cells the two agree, so the usual examples are unaffected; with elongated cells, e.g. a drawn geometry
extruded along y with a different cell size there, or `slump_problem` with unequal `L ./ nel`, the step
is set by the largest cell and may exceed the stable step of the smallest. Found while writing
`estimate_dt`, which mirrors `get_dt` as it is (and is tested equal to it), so fixing
`get_dt` needs the same change there.
