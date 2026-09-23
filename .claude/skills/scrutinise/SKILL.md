---
name: scrutinise
description: Scrutinise newly added or changed code on the current branch against main. Checks new types, methods, changed signatures, and helpers for necessity, correctness, clarity, consistency, robustness, and minimality. Reviews new tests for gap coverage, overlap with existing tests, and minimality. Invoke with /scrutinise.
tools: Bash, Glob, Grep, Read, Edit, Write
---

# Scrutinise

Review newly added or changed code on the current branch against `main`.

## 1. Gather the diff

```bash
git diff main...HEAD --stat
git diff main...HEAD -- src/ test/
```

Read every changed file in full, not just the hunks.

## 2. Source review (`src/`)

The package layout is: `src/boot/` (includer + hand-ordered `types/`),
`src/home/core/common/` (shared ignite/shpfun/topology kernels), `src/home/core/solver/
explicit/` and `.../dynamic_relaxation/` (the two solver paths), `src/home/init/`
(setup_mesh/setup_mpts/setup_problem/setup_basis), `src/home/api/solver/` (get_default,
get_solver, cli). Full detail: `.claude/docs/architecture.md`.

For each new/changed type, method, or helper check:

- **Necessity** — does existing infrastructure already provide this? (e.g. `Basis`
  already owns all connectivity — a new type shouldn't duplicate `e2n`/`p2n`; a new
  dispatch axis should extend the existing "config string → marker instance →
  multiple dispatch" pattern rather than adding a parallel branch.)
- **Correctness** — type parameters consistent with the surrounding struct
  (`Point{T1,T2,D,CM,...,ST,SC,SK}`, `Basis{T1,T2,D,NN,K,TR}`); `NN`/connectivity stays
  on `Basis`, never `Point`/`Mesh` (see `conventions.md`); a stored tensor
  (`mpts.s.σᵢⱼ`/`ϵᵢⱼ`) is read via `get_voigt`/`get_tensor`, never `.dev`/`.p` directly;
  `basis.trsfr`/`basis.which` stay orthogonal — a change to one shouldn't implicitly
  constrain the other.
- **GPU/KernelAbstractions compatibility** — anything running inside a `@kernel` body
  avoids scalar indexing into device arrays, dynamic dispatch, or CPU-only control
  flow; `@index(Global)` is wrapped `T1(p)` before use; bare float literals are wrapped
  `T2(...)` in any function generic over the float type (see `gotchas.md`).
- **Precision (32-bit) correctness** — if the change touches anything generic over
  `T1`/`T2`, has it actually been exercised under `dtype=(;T0=(Int32,Float32),...)`, or
  only under the default `(Int64,Float64)` where these bugs silently pass?
- **Clarity & consistency** — unicode-subscript tensor-rank convention followed
  (`ᵢⱼ` for a full 2nd-order tensor, `ᵢ` for a genuine vector — see `conventions.md`);
  naming and idiom match the surrounding file.
- **Minimality** — no dead code, no speculative generality, no new config knob or
  abstraction sized for a hypothetical future need (see the repo's "less is better"
  motto) — prefer running `/minimise` on the diff before or alongside this review.

## 3. Test review (`test/`)

- Does each new test cover a genuine gap, or does an existing
  `test/testset/test_*.jl` file already exercise it?
- Are new tests in the right file (`test_basis.jl`, `test_column.jl`,
  `test_performance.jl`, `test_workflow.jl`) rather than a new file?
- Is a refactor verified the way `conventions.md` requires — an actual
  `slump_problem` → `elastoplasm!(...; workflows=[...])` run with real post-run field
  values inspected, not just a clean-load/`success=true` check?
- If the change touches `test_performance.jl`-benchmarked kernels, was it checked
  against `test/dataset/performance_baseline.jld2`, and is any intentional shift
  documented (see `.claude/docs/conventions.md`)?

## 4. Output

Report findings grouped by file, each tagged with the failed criterion (necessity /
correctness / GPU compatibility / precision / clarity / minimality). Propose concrete
simplifications. Apply fixes only when asked.
