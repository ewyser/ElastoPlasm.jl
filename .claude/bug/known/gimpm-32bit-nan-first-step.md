# 32-bit GIMP run fails at the first step with `InexactError: trunc(Int32, NaN)`

**Status: open, root cause not investigated.** Workaround: run GIMP in 64-bit (the default).

Under `dtype = (; T0=(Int32, Float32), bits=Int32(32), precision="32-bit precision (or single precision)")`,
the stock 2D slump (`slump_problem([64.1584, 16.0396], [40, 10]; cli()..., dtype=...)`) with
`basis = (; get_default().basis..., which="gimpm")` throws `InexactError: trunc(Int32, NaN)` in
`elastodynamic!`'s first time step. A NaN reaches an `Int32` conversion, presumably an element
index computed from a particle position or domain (`p2e`/`tplgy!`). Pre-existing: the same call
fails before per-region elastic constants were added.

In 32-bit, the same slump runs with the default basis, with `material.plastic="VM"`, with
`stab.locking=false`, and with `bcs.dirichlet=[:roller :roller; :fixed :free]`, so the failure is
specific to GIMP. It was found while testing the `dtype.precision` option.
See `gotchas.md` for the two 32-bit patterns already fixed (`@index(Global)` and bare literals).
