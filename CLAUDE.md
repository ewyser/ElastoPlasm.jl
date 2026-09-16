# ElastoPlasm.jl — Notes for Claude

Working notes on this repo's architecture and conventions, kept up to date as the
codebase evolves. Read this first when picking up work here.

This file is the overview: purpose, structure, core method, and git conventions. Detail
lives in `.claude/docs/` (auto-loaded below via `@` imports, so nothing here is lost by
splitting it out) and in `bug/` (currently-open bugs under `bug/known/`, resolved bugs
with full investigation history under `bug/fixed/` — one file per issue, referenced on
demand rather than imported wholesale).

## Project shape

MPM/geomechanics solver in Julia. `src/home/` is loaded via a custom recursive includer
(`src/boot/include.jl`, `collect_and_include_jls`) that walks the directory tree and
includes every `.jl` file — files in a directory load before its subdirectories, both
sorted alphabetically, so load order is deterministic there. A failed include only
`@warn`s, it doesn't crash — so a stale/broken file can silently go unused; always
check for `@warn "Failed to include..."` after `using ElastoPlasm` when touching files
under `home/`.

`src/boot/needs/types/` is loaded differently: `boot.jl` includes its files through an
explicit, hand-ordered list rather than `superInc`/`collect_and_include_jls`, because
alphabetical order doesn't match dependency order there (e.g. `tensor.jl`'s
`AbstractStrain`/`AbstractStress` must exist before `problem/lagrangian.jl`'s `Point`
constrains against them). Each abstract type is defined directly in the file that
implements its concrete subtype(s) (no central `abstract.jl`) — see
`.claude/docs/architecture.md` for the directory layout. Adding a new file under
`types/` means adding it to that explicit list in `boot.jl`, in the right dependency
position; it is NOT picked up automatically the way `home/` files are.

## Core method, at a glance

- Construction order is always **mesh → material points → basis** (basis depends on
  both). Function argument order is always **`(mpts, mesh, basis, ...)`**.
- Three core types: `Mesh` (Eulerian grid), `Point` (Lagrangian material points),
  `Basis` (owns all connectivity + shape-function kind + transfer-scheme dispatch).
  Full field-by-field detail, transfer-scheme dispatch, and the explicit vs.
  dynamic_relaxation solver split: `.claude/docs/architecture.md`.
- Config/behavior knobs (`basis.which`, `strain.deform`, `stab.locking`, `plast.*`,
  `nonloc.*`, ...) live on one `NamedTuple` built by `get_default()`/`get_solver` and
  merged shallowly per-section. Full reference and the standard end-to-end run:
  `.claude/docs/operations.md`.
- Persistence layout (`ic["problem"]`/`ic["basis"]`/`cfg["solver"]` in JLD2):
  `.claude/docs/persistence.md`.
- 32-bit precision and StaticArrays allocation-elision gotchas:
  `.claude/docs/gotchas.md`.
- House rules (unicode subscript tensor-rank convention, refactor-verification
  discipline, `dump/` flushing on load, etc.): `.claude/docs/conventions.md`.
- Known-follow-up design work (not bugs): `.claude/docs/planned-improvements.md`.
- Open bugs: `bug/known/*.md`. Resolved bugs, kept for archaeology: `bug/fixed/*.md`.
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
`bug/known/`/`bug/fixed/` against what actually happened, and update whichever files
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
