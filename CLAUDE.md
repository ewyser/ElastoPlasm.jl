# ElastoPlasm.jl — Notes for Claude

Working notes on this repo's architecture and conventions, kept up to date as the
codebase evolves. Read this first when picking up work here.

This file is the overview: purpose, structure, core method, and git conventions. Detail
lives in `.claude/docs/` (auto-loaded below via `@` imports, so nothing here is lost by
splitting it out) and in `.claude/bug/` (currently-open bugs under `.claude/bug/known/`, resolved bugs
with full investigation history under `.claude/bug/fixed/` — one file per issue, referenced on
demand rather than imported wholesale).

## Project shape

MPM/geomechanics solver in Julia. `src/home/` is loaded via a custom recursive includer
(`src/boot/include.jl`, `collect_and_include_jls`) that walks the directory tree and
includes every `.jl` file — files in a directory load before its subdirectories, both
sorted alphabetically, so load order is deterministic there. A failed include only
`@warn`s, it doesn't crash — so a stale/broken file can silently go unused; always
check for `@warn "Failed to include..."` after `using ElastoPlasm` when touching files
under `home/`. The warning appears only while the package **precompiles** (the first
`using` after an edit): later sessions load the cached image, which silently lacks the
broken file's definitions (an `UndefVarError` for something that exists in the source).
To see it again, `include` the file directly or force recompilation. Adding or removing a
file recompiles (`superDir` registers every directory with `include_dependency`; before
that, a new file was silently missing from the cached image).

`src/boot/needs/types/` is loaded differently: `boot.jl` includes its files through an
explicit, hand-ordered list rather than `superInc`/`collect_and_include_jls`, because
alphabetical order doesn't match dependency order there (e.g. `tensor.jl`'s
`AbstractStrain`/`AbstractStress` must exist before `problem/lagrangian.jl`'s `Point`
constrains against them). Each abstract type is defined directly in the file that
implements its concrete subtype(s) (no central `abstract.jl`) — see
`.claude/docs/architecture.md` for the directory layout. Adding a new file under
`types/` means adding it to that explicit list in `boot.jl`, in the right dependency
position; it is NOT picked up automatically the way `home/` files are.

`lib/EasyDraw/` is EasyDraw.jl (https://github.com/ewyser/EasyDraw.jl), a package and repository of its own,
included here as a **git submodule**, usable without ElastoPlasm (a friend's finite-difference
code, finite elements): the shapes and `Zone`, layouts, the drawing window (QML, a dependency loaded by the first window),
sampling and the partition for Gmsh. A change to EasyDraw is committed and pushed in `lib/EasyDraw` (its repository),
then ElastoPlasm commits the new submodule commit (`git add lib/EasyDraw`); see "Working with the submodules"
in `.claude/docs/gotchas.md`. ElastoPlasm depends on it by path (`Pkg.develop`, recorded in
the Manifest; Julia 1.10 has no `[sources]`) and is one host of its window (`get_catalog`, `case_hooks`); a
`Case` is built from an EasyDraw `Layout` and its `Settings` (`Case(layout, settings)`), so a case file is a layout
file with settings and a `Case` works wherever EasyDraw takes a layout, with the same points to the bit; its tests run from `test/testset/test_easydraw.jl`. Changing
EasyDraw's dependencies needs `Pkg.resolve()` in ElastoPlasm's environment. See "Drawing a case" in
`.claude/docs/operations.md`.

`lib/UnifiedBackend/` is UnifiedBackend.jl (https://github.com/ewyser/UnifiedBackend.jl), also a submodule developed by
path: the registry of where kernels run, `Platform`s with their KernelAbstractions `backend` and `arraytype` (the host,
any CPU, registered at load; GPUs registered by its CUDA/ROCm extensions when CUDA.jl/AMDGPU.jl loads and works).
`build_solver` stores `UnifiedBackend.select_platform(backend.select)` (`"host"`/`"device"`, the first GPU or else the
host) as `solver.backend.exec`, activated (`activate!`); `exec_backend(solver)` is its backend. ElastoPlasm does
`import UnifiedBackend` and qualifies its names. It replaced ElastoPlasm's own, older copy of the same code
(`src/boot/needs/backend.jl`, `ext/`); version 0.1.0 replaced the earlier API (`Host`/`Device`, `get_backend`,
`select_execution_backend`, menus).

## Core method, at a glance

- A problem starts from a `Case`: an EasyDraw `Layout` (the geometry: domain, grid, named zones) and its `Settings`
  (the values: zones' constants by name, material defaults, sides, solver options, run settings), plain data drawn or
  written, that `Case(layout, settings)` checks and reads as typed values; set up by `setup_simulation` and stored as one
  TOML (`save_case`/`load_case`; geometry and settings may be separate files); every example is a case constructor
  (`slump_case`, …), and the drawing window (EasyDraw's) edits both. See "Describing a case" in
  `.claude/docs/operations.md`.
- Construction order is always **mesh → material points → basis** (basis depends on
  both). Function argument order is always **`(mpts, mesh, basis, ...)`**.
- Three core types: `Mesh` (Eulerian grid), `Point` (Lagrangian material points),
  `Basis` (owns all connectivity + shape-function kind + transfer-scheme dispatch). What is solved, a `Case`
  (geometry with its physics), and how, a `Solver` (`Solver(case; overrides…)` from a case's settings, or
  `Solver(; dim, …)`: numerics, run settings, kernels; it keeps nothing of a case), meet in
  `Problem(case, solver)` (the three, with `Time` and the solver); a `Simulation` is a problem ready to run in a
  folder (the case's), run by `execute!`. One solver solves several cases of its dimension.
  The type diagram, field-by-field detail, transfer-scheme dispatch, and the explicit vs.
  dynamic_relaxation solver split: `.claude/docs/architecture.md`.
- Config/behavior knobs (`basis.which`, `material.elastic`, `material.plastic`, `stab.locking`,
  `nonloc.*`, ...) live on one `NamedTuple` built by `get_default()` and
  merged per field (`merge_config`); the physics sections (`material`, `grf`, `bcs`) are a case's, the rest the
  solver's. Full reference and the standard end-to-end run:
  `.claude/docs/operations.md`.
- Persistence layout (one `"simulation"` key, the whole `Simulation`, in JLD2):
  `.claude/docs/persistence.md`.
- GPUs: the explicit path runs on a device (`backend = (; select = "device")`, a GPU UnifiedBackend registered): the
  containers' array types are type parameters, the kernels are built on the solver's backend, the workflows move
  the problem to the device and back (`core/device.jl`). Tested on JLArrays (`test_device.jl`), validated on an
  NVIDIA machine by `test/cuda/validate_cuda.jl`; dynamic relaxation and thermal still host only. See "Running on a
  GPU" in `.claude/docs/operations.md`.
- 32-bit precision and StaticArrays allocation-elision gotchas:
  `.claude/docs/gotchas.md`.
- House rules (unicode subscript tensor-rank convention, refactor-verification
  discipline, `dump/` emptied only by `clean_dump()`, etc.): `.claude/docs/conventions.md`.
- Known-follow-up design work (not bugs): `.claude/docs/planned-improvements.md`.
- Open bugs: `.claude/bug/known/*.md`. Resolved bugs, kept for archaeology: `.claude/bug/fixed/*.md`.
  **When picking up an open bug: create a new branch off the current branch, named so
  the fix is identifiable (e.g. `fix-volumetric-locking-zero-mass-guard`), rather than
  fixing in place on an unrelated branch.**
- Specialized subagents (`.claude/agents/`, invoke explicitly via the Agent tool):
  `agent-user` (usability/DX audit), `agent-mpm-specialist` (MPM/geomechanics science
  correctness against published literature in `../refs/ElastoPlasm/`),
  `agent-computer-science-specialist` (Julia dispatch/KernelAbstractions/StaticArrays/
  CUDA performance). Each is read-only — reports findings rather than editing code.
  For a decision spanning more than one of these lenses, invoke `agent-arbiter`
  instead — it dispatches the relevant specialists in parallel and synthesizes their
  findings into one tradeoff-aware recommendation, rather than you reconciling three
  separate reports yourself.

@.claude/docs/architecture.md
@.claude/docs/operations.md
@.claude/docs/persistence.md
@.claude/docs/gotchas.md
@.claude/docs/conventions.md
@.claude/docs/planned-improvements.md

## End of every task

Before considering a task done, review the six `.claude/docs/*.md` files and
`.claude/bug/known/`/`.claude/bug/fixed/` against what actually happened, and update whichever files
have gone stale — a fixed bug, a newly found bug, a design decision, or identified
follow-up work. This is a standing expectation for every task, not something that only
happens when `/log-session` is invoked by name; use that skill when the update is
substantial enough to want its file-matching/status-line discipline, but don't treat
its absence as license to skip the review on a smaller change.

## Commit message conventions

This repo follows Conventional Commits: `<type>: <Summary>`, summary capitalized, no
trailing period. Pick the type by what actually changed, not by what the change was
*for*:

- `feat:` — new functionality or capability that wasn't there before.
- `fix:` — corrects behavior that was wrong (a bug, a crash, an incorrect result).
- `refactor:` — restructures existing code without changing external behavior (e.g.
  renaming, moving logic between files, simplifying a dispatch pattern).
- `chore:` — maintenance with no production/behavior impact: build tooling, dependency
  bumps, config, repo housekeeping, formatting.
- `docs:` — documentation-only changes (README, CLAUDE.md, docstrings).
- `test:` — adding or updating tests only, no source changes.
- `perf:` — a performance improvement, not a correctness fix.

`git log --oneline` in this repo also has some older/looser types in its history
(`misc:`, `bug:`) — treat those as historical, not part of the current convention;
prefer the list above for new commits.
