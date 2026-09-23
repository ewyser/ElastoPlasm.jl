## Conventions / house rules

- **Unicode subscripts encode tensor rank**: a full 2nd-order tensor field/variable
  (stress, strain, deformation gradient, ...) uses the double-index `ᵢⱼ` suffix (e.g.
  `σᵢⱼ`, `τᵢⱼ`, `ϵᵢⱼ`, `Fᵢⱼ`, `Bᵢⱼ`, `Dᵢⱼ`); a genuine vector quantity (velocity,
  displacement component, a Voigt-notation vector like `get_voigt(...)`'s result)
  uses the single-index `ᵢ` suffix. Never use `ᵢ` on a field that actually stores a
  full tensor, or `ᵢⱼ` on a plain Voigt vector — both misread as the wrong rank.
- `NN` (nodes-per-element count) belongs to `Basis`, never to `Point` or `Mesh` — it's
  topology, not point/mesh state.
- Prefer `SVector`/`SMatrix` (StaticArrays.jl) over runtime-sized arrays for per-point
  fixed-size data.
- When converting one solver path to a new convention, apply the identical mechanical
  pattern already validated on the sibling path (explicit ↔ dynamic_relaxation) rather
  than reinventing it — these two solver directories are meant to stay structurally
  parallel.
- Verify refactors by actually running the `slump_problem` → `elastoplasm!(...;
  workflows=[...])` pipeline end-to-end and inspecting real post-run field values, not
  just a clean-load check or `success=true` — a file can `include` cleanly, a workflow
  can return `success=true`, and the run can still be functionally broken (wrong
  kwarg name, stale dict key, matrix-style indexing into what's now a `Vector{SVector}`,
  or — as happened once — silently operating on the untouched pre-simulation state
  because `elastoplasm` was used instead of `elastoplasm!`).
- `test/testset/test_performance.jl` benchmarks core kernels directly via
  `solver.cairn.*`; keep it in sync with the `Cairn` dispatch-table shape whenever that
  shape changes.
- `using ElastoPlasm` flushes (deletes) everything under `dump/` on load if that
  directory already has contents. Never run more than one `julia --project=. -e '...'`
  invocation against this repo concurrently while any of them is still writing to
  `dump/` — a second session's `using ElastoPlasm` will delete the first session's
  still-in-progress output out from under it.
- When doing a mechanical rename across a codebase (e.g. via sed/substring match),
  verify each occurrence semantically before trusting it — a blind rename can catch
  things that only textually match but mean something different (see the `σᵢ`/`τᵢ`
  entry in `planned-improvements.md` for a concrete example: local Voigt-vector
  variables getting swept up in a tensor-field rename).
