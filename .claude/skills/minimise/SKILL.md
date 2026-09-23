---
name: minimise
description: Prune a bug fix or new tests down to the smallest correct diff through multiple elimination passes. Use before committing any fix or test addition, in the spirit of this repo's "less is better — no AI slop" motto.
---

# Minimise

The goal is to remove every line that is not strictly required for correctness, then
verify the result still passes the relevant tests. This operationalizes the repo's own
standing motto (`.claude/docs` header / pinned memory: "less is better — no AI slop")
as a repeatable pass rather than a value someone has to remember to apply by eye.

## Process

Repeat the following until no further reductions are possible:

1. **Read the diff.** Run `git diff HEAD` (or `git diff --cached` if staged) and read
   every changed file in full.

2. **Challenge each change.** For every changed line ask:
   - Would removing this line cause a test to fail or a bug to reappear?
   - Is this a cleanup, rename, refactor, or comment that is not load-bearing?
   - For new tests: does an existing `test/testset/test_*.jl` file already cover this
     behaviour?  If so, drop the new test entirely.
   - Is this a new abstraction (type, config knob, dispatch layer) solving a narrower
     problem than its surface area suggests? Prefer the narrower, more direct fix.

3. **Remove non-essential changes.** Delete anything that does not answer "yes" to the
   first question above. Prefer shrinking an existing case over adding a new one.

4. **Run the focused tests.** `test/runtests.jl` interactively selects which
   `test/testset/test_*.jl` files to run (auto-runs all under `GITHUB_ACTIONS=true`);
   it swallows exceptions inside its own `try/catch`, so `include()` the relevant file
   directly when you need the real stacktrace — see `.claude/docs/operations.md` for
   the exact preamble. Run the full suite
   (`julia --project=. -e 'using Pkg; Pkg.test()'`) before the final report.

5. **Repeat** from step 1 until a full pass produces no further removals.

## Heuristics

- A one-line fix is better than a five-line fix.
- A new test case added to an existing `@testset` is better than a new test file.
- Comments and blank lines added alongside a fix are not load-bearing; remove them
  unless they explain something non-obvious (a hidden constraint, a workaround, a
  subtle invariant — see `.claude/docs/conventions.md`).
- Helper functions introduced solely for the fix are a red flag; inline them unless
  the surrounding file already uses that granularity of decomposition.
- Changes under `src/home/script/example/` are demos, not library code: keep them out
  of a core-kernel bug-fix diff unless the example itself is the regression check.
- A refactor-driven change that touches both `src/home/core/solver/explicit/` and
  `src/home/core/solver/dynamic_relaxation/` should apply the identical mechanical
  pattern to both, not a bespoke variant per path (see `conventions.md`).

## When to stop

Stop when every remaining line answers "yes" to: *if I remove this, the targeted bug
reappears or the targeted test fails*. At that point report the final diff and suggest
committing — but do not commit without being asked (see `CLAUDE.md`'s git conventions).
