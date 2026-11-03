## Conventions / house rules

- **Unicode subscripts reflect how the data is stored**: `ᵢ` for anything stored as a vector —
  including stress and strain, stored as a scalar plus a 3D Voigt deviator (`σᵢ`, `τᵢ`, `ϵᵢ`) — and
  `ᵢⱼ` for anything stored as a matrix (`Fᵢⱼ`, `ΔFᵢⱼ`, `∇vᵢⱼ`, `Bᵢⱼ`, `Dᵢⱼ`, or a local
  `σᵢⱼ = get_tensor(...)`). This replaced an earlier rule tied to physical tensor rank. Don't
  give a temporary, differently-shaped object a stored field's name (e.g. a 3×3 embedding of
  `ΔFᵢⱼ`); convert visibly inside the function, under a comment saying what it models.
- **Writing a type** (one way, asked for by the user; the pass is in `planned-improvements.md`, "Type writing"):
  - `T1` is always the index type and `T2` the float type, in type parameters and in methods' `where` clauses (no
    plain `T` for either); `D` the dimension.
  - Bound only the parameters code dispatches on as kinds (`CM<:AbstractConstitutiveModel`, `ST`, `EL`, `K`, `TR`,
    GIMP's domain); leave number, array and part types unbounded (bounds referring
    to other parameters made compilation three times slower, `architecture.md`). Array and part parameters come after
    the dispatched ones, in field order, with a constructor taking them from the arguments and checking what must agree
    (`Mesh{T1,T2,D}(…)`; `Problem(case, solver)` builds its parts in the solver's `T1`, `T2`, `D`; `Point(…)` takes
    its dispatched parameters from its solid phase).
  - Parameter lists without spaces (`{T1,T2,D}`); fields as `name ::T`, aligned.
  - An abstract type only when code dispatches on it, named after the role of the field it types
    (`AbstractBasisKind` for `Basis`'s `kind`, `AbstractTransfer`, `AbstractDomain`); its concrete markers named
    `<Name><Role>` (`StdTransfer`, `StretchDomain`). A marker type only when something dispatches on it: the solver's
    `solution` stays text with one check (`explicit_only`), since nothing dispatched on its markers (the user's
    decision). Where an option's values name types needed both ways (name → type to build, type → name to report), one
    function returns the table of both, a NamedTuple keyed by the option's values in their listed order
    (`get_constitutive_laws()`, read by `table_value`, `laws(mpts)` the other way; no per-name lookup function); a
    NamedTuple, not a `Dict`, for the order and like the other in-memory tables (`Dict`s are the file form).
    `get_option`'s lists derive from it. `get_basis`, `get_transfer` and `get_domain` build their markers from a name
    with arguments, one way only, as branching functions.
- **Keep things independent of each other wherever possible** (the user's standing goal): a part knows only what it
  needs, and the direction of what may know what is explicit. EasyDraw knows nothing of ElastoPlasm; a `Layout` knows
  no values; `Settings` hold values by zone name with no meaning given to them; a `Case` (EasyDraw's) checks only that
  the two agree; ElastoPlasm gives the values their meaning (`case_physics`, `Solver(case)`); a `Solver` holds nothing
  of a case. In EasyDraw's window, `Geometry` knows no values and `Values` know zones only by id; a page knows the
  parts it edits or shows; the frame knows no page; which pages a window has is the opening function's choice. A new
  page is a type, its methods and a QML file in `src/pages/<page>/`, nothing elsewhere.
- **Every table is a function, never a `const`** (agreed with the user, like `get_default()`, `get_option()`):
  `get_<what>()`, plural for a collection (`get_workflows()`, `get_boundary_kinds()`, `get_zone_properties()`,
  `get_physics_sections()`, `get_constitutive_laws()`), singular for one record (`get_material()`, `get_run()`). A table
  read in one function only is a local there (the catalog's page order). `const` stays for what is not a table: the
  source path (`SRC`); in EasyDraw, the open window's editor (`EDITOR`, a `Ref`) and its QML file. EasyDraw follows
  the rule in its own naming (nouns: `side_names(D)`, `canvas_padding()`, `qml_calls()`).
  - Every type's docstring in the house layout: the signature with its parameters, a summary, `# Fields` (each field's
    type, values and units; then what the type parameters are), `# Constructors`, `# Throws`, `# Examples`,
    `See also`, only the sections that apply. A field nothing reads says so.
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
- Verify refactors by actually running the `setup_simulation(slump_case(...))` → `execute!(jld2)`
  pipeline end-to-end and inspecting real post-run field values (`load_simulation(jld2).problem`),
  not just a clean load — a file can `include` cleanly and a run can finish and still be
  functionally broken (wrong kwarg name, stale dict key, matrix-style indexing into what's now a
  `Vector{SVector}`, or — as happened once with the old `elastoplasm(jld2)`, which discarded its result — silently
  operating on the untouched pre-simulation state; `execute(jld2)` returns its result but leaves the
  file, so read a run's results from what it returns, or run `execute!`). For plumbing changes, compare final states bit for bit
  against the parent commit, single-threaded (`gotchas.md`).
- `test/testset/test_performance.jl` benchmarks core kernels directly via
  `solver.cairn.*`; keep it in sync with the `Cairn` dispatch-table shape whenever that
  shape changes.
- `dump/` is never emptied automatically; `clean_dump()` empties it on request. (Until
  `merge-hdf5-output` it was emptied by top-level code in `boot.jl`, i.e. whenever the package
  precompiled: the first `using` after a source change deleted every output, also those of another
  session still running; the user chose no automatic emptying over a prompt in `__init__`, which
  would hang scripts, tests and CI.) `boot.jl` only creates it, since git does not track it.
  `setup_simulation` overwrites the files of its own simulation and removes its stale output, so
  older simulations' folders stay until `clean_dump()`.
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
