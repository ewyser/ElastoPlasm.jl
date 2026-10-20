## Conventions / house rules

- **Unicode subscripts reflect how the data is stored**: `ᵢ` for anything stored as a vector —
  including stress and strain, stored as a scalar plus a 3D Voigt deviator (`σᵢ`, `τᵢ`, `ϵᵢ`) — and
  `ᵢⱼ` for anything stored as a matrix (`Fᵢⱼ`, `ΔFᵢⱼ`, `∇vᵢⱼ`, `Bᵢⱼ`, `Dᵢⱼ`, or a local
  `σᵢⱼ = get_tensor(...)`). This replaced an earlier rule tied to physical tensor rank. Don't
  give a temporary, differently-shaped object a stored field's name (e.g. a 3×3 embedding of
  `ΔFᵢⱼ`); convert visibly inside the function, under a comment saying what it models.
- `NN` (nodes a point's shape functions reach) belongs to `Basis`, never to `Point` or `Mesh` — it's
  topology, not point/mesh state. Read it as `stencil_size(basis)` (the `NN` type parameter, so
  `for nn ∈ 1:stencil_size(basis)` has a compile-time bound). `MeshProperties.nn` and `Point.nn`,
  copies of it, were removed; the switch was bit-identical and speed-neutral (2 s slump, bsmpm and
  gimpm, single thread: 1.10 s before and after).
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
- `dump/` is emptied when the package **precompiles** (`rootflush` is top-level code in
  `boot.jl`, so it runs while the cache is built, not at every `using`): the first `using
  ElastoPlasm` after a source change deletes everything under `dump/`; later sessions loading the
  cached image leave it alone (checked: a probe file survived a plain `using`). So outputs survive
  between sessions until the next edit, and stale files (e.g. a previous run's `_output.h5`) can
  still be there; `setup_simulation` removes the output of the problem it rewrites. Never run more
  than one session against this repo while one of them may recompile and another is writing to
  `dump/`.
- When doing a mechanical rename across a codebase (e.g. via sed/substring match),
  verify each occurrence semantically before trusting it — a blind rename can catch
  things that only textually match but mean something different (see the `σᵢ`/`τᵢ`
  entry in `planned-improvements.md` for a concrete example: local Voigt-vector
  variables getting swept up in a tensor-field rename).
- **Stress and strain are full 3D tensors in every dimension, stored split.** Stress is `p`
  (pressure, positive in compression) plus a trace-free Voigt `dev`; strain is `vol` (trace) plus
  a trace-free Voigt `dev` (engineering shear). Constitutive code (elastic laws, return maps) works
  on `.p`/`.vol`/`.dev` directly and must pass a trace-free `dev` when it builds one. Kinematics
  meets the split only through the 3×3 constructors (tensor → split) and `get_tensor`
  (split → tensor). Write constitutive code once, in 3D: kinematics stays D×D and enters
  through one carry constructor, `S(M, t)` (new in-plane D×D block `M`, out-of-plane components
  carried from `t`, for any stress or strain type), and never add a 2D-specific constitutive method.
