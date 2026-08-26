# Ongoing task state — typed strain/stress tensor work & follow-ups

Local, gitignored (as of this session — `.claude/` is now in the repo's own
`.gitignore`, not just the user's global git ignore, which previously only covered
`.claude/settings.local.json` specifically) trace of in-progress Claude Code work on
this repo. Check this file, and `.claude/plans.md` (queued-but-not-started plans), at
the start of a session here.

**Note on this file's history**: an earlier, more detailed version of this file
existed earlier in this same session (tracking the typed-tensor port round by round)
but went missing partway through — found gone, cause unknown, when checked again
later in the session. Given `.claude/` was never actually gitignored until this
rewrite (see above), it's possible something touched it; no direct evidence found.
Rebuilt from conversation context below rather than left blank.

## Done, merged to `dev`

The full typed strain/stress tensor port and its follow-up cleanups, on branch
`feat-typed-strain-stress-tensors` (5 commits: `a026b4f` tensor port,
`fbc64ec` Δλ fixes, `dfc2811` `AbstractElasticity` removal, `67e5297` `get_voigt`
convention fix + `mutate`/`_mutate` retirement, `941e8d4` constructor/comment
cleanup) — merged into `mimic-implementation-logic-of-ample`, then into `dev`
(merge commits `ff6134a`, `78838bf`), and `mimic` fast-forwarded to match `dev`. One
follow-up doc commit landed directly on `mimic`/`dev` after the merge: `2d474c5`
(`docs: Confirm test_workflow.jl baseline holds after tensor-port merge to dev`).

All of this is independently verified (not just self-reported) at multiple points
this session: bit-exact/machine-epsilon numerical equivalence checks against pre-change
math, zero-allocation checks, the full `test_workflow.jl` 96-case sweep unchanged at
71 passed/25 failed throughout, `test_collapse.jl` 4/4, `test_basis.jl` passing,
`test_performance.jl` showing real improvement (-17.2% memory/-5.4% allocs/0.0% time
vs. the `dev`-generated baseline) — see CLAUDE.md's "Typed strain/stress tensor
storage" and related entries for the full detail, since those are committed and won't
go missing the way this file did.

`fancy-implementation` (the source branch these ideas were ported from) still exists
locally and on `origin` — user was weighing deleting it since everything useful from
it is now ported; no decision recorded on whether/when to actually delete.

## Done, merged to `dev` (update — supersedes the "Queued" section below)

The transfer-scheme-dispatch plan (`.claude/plans.md`) is **done**, not queued:
`Basis` gained a `transfer::TR<:AbstractTransfer` field (`StdTransfer`/`TpicTransfer`/
`ApicTransfer{T,D,L}`, the latter carrying `Bᵢⱼ`/`Dᵢⱼ`, moved off `Point` — `Δnp` was
dead code, deleted). `p2n!`/`Bij` now dispatch on `TR` instead of a `Cairn`
string-branch. Landed as `feat-dispatch-transfer-scheme-on-basis`, then a follow-up
config restructuring (`refactor-fold-transfer-config-into-basis-stab`: folded
`transfer.trsfr`/`.C_pf` into `basis`, `transfer.musl` into `stab`, removed the
`transfer` config section entirely) — both merged into `mimic-implementation-logic-
of-ample` then `dev`, `mimic` fast-forwarded to match. See CLAUDE.md's "Transfer
scheme dispatch via `Basis`'s `TR` type parameter" entry for full detail.

Also done, same merge lineage: DP/J2 `retmap` kernel unification (dispatched on
`Point`'s `CM`/`ST`, new `VonMises` constitutive model), and a `shpfun!`/`MLSBasis`
dispatch-ambiguity fix (`test_workflow.jl`'s 96-case sweep baseline is 71 passed/25
failed, confirmed unchanged across all of the above) — see CLAUDE.md Known bugs.

**2026-08-26**: did a full compression pass on CLAUDE.md (1373→~750 lines, commit
`f964427` on `mimic-implementation-logic-of-ample`) — collapsed multi-layer
"done, then re-verified" narratives for now-stable/merged features into current-state
descriptions, kept every load-bearing gotcha and the Known-bugs section's diagnostic
content intact, fixed a few stale `transfer.trsfr` references left over from the
config restructuring above.

## Done: performance comparison vs. `MaterialPointSolver.jl`

**2026-08-26, completed.** Ran a matched CPU comparison against
`LandslideSIM/Archive_MaterialPointSolver.jl_paper`'s `02_soil_collapse/2d_collapse.jl`
(the vendored paper-version solver, not current `main`/v0.5.1 — see below for why),
`device=:CUDA`→`:CPU`. Results (also given to the user in chat):

| | Particles | Iterations | Wall time | Throughput |
|---|---|---|---|---|
| MaterialPointSolver.jl (uGIMP, DP, MUSL, `:CPU`) | 12,800 | 71,215 | 458.0 s | 0.50 µs/particle/iter |
| ElastoPlasm `elastodynamic!` (bsmpm, `std`, no plasticity) | 10,133 | 1,137 | 12.0 s | 1.04 µs/particle/iter |
| ElastoPlasm `elastoplastic!` (bsmpm, DP + F-bar + non-local reg.) | 10,133 | 640 | ~193.1 s | 29.8 µs/particle/iter |

Conclusion given to and confirmed by the user: ElastoPlasm's raw explicit-dynamics
kernel is competitive (~2x slower per particle-step, not an order of magnitude), and
the much larger `elastoplastic!` gap isn't a fair like-for-like — that phase runs
F-bar locking correction + non-local regularization, both default-on in ElastoPlasm
and absent from the comparison run, not pure overhead.

**Follow-up, same day**: re-ran with `stab.locking=false`+`nonloc.status=false` (user
request). Result: `elastoplastic!` dropped from 193.1s→**5.0s** for essentially the
same iteration count (640→637), i.e. **0.78 µs/particle/iter — only ~1.6x slower**
than MaterialPointSolver.jl's combined MUSL+DP step (0.50 µs/particle/iter), not the
~60x the locking+non-local-inclusive number showed. Confirms those two features
account for nearly the entire earlier gap, not the DP retmap/basis evaluation itself.

Also generated a visual (before/after scatter plot of a small `2d_collapse.jl` run,
3200 particles, `:CPU`) to show what the reference solver's own output looks like —
saved to the session scratchpad only (`mps_bench/mps_collapse_result.png`), not
copied anywhere durable (user didn't ask to keep it).

**Important scope clarification surfaced by the user, worth remembering if this is
revisited**: the comparison only ever matched *particle count/order of magnitude and
material config* between the two sides — it did **not** match initial geometry.
MaterialPointSolver.jl's `2d_collapse.jl` releases a rectangular block (a "collapse");
ElastoPlasm's `slump_problem` starts from a pre-carved sloped/wedge shape (a "slump").
The plotted collapse image is `2d_collapse.jl`'s own shape only — it was never
generated for the ElastoPlasm run. A true same-geometry comparison (offered to the
user, not yet requested) would need either ElastoPlasm's `collapse_problem`
(rectangular block, matches MaterialPointSolver.jl's own geometry) in place of
`slump_problem`, or a hand-built wedge IC on the MaterialPointSolver.jl side.

**User's final assessment (2026-08-26), worth keeping verbatim as the session's
conclusion on this task**: "ElastoPlasm holds, especially given the fact that it is
more generic, there is no high-performance oriented implementation logic" — i.e. the
~1.6-2x gap (once locking/non-local are excluded from the comparison) is considered a
good result given MaterialPointSolver.jl is purpose-built for memory-throughput
performance (per its own paper) while ElastoPlasm prioritizes generality (pluggable
basis/transfer/constitutive-model/strain-type dispatch) over raw-throughput tuning.

## Done: real granular-collapse example added; misnamed old one renamed

**2026-08-26, same day, follow-up.** User's next request: implement an actual granular
block-collapse example (the comparison above used `slump_problem`'s pre-sloped geometry
on the ElastoPlasm side vs. MaterialPointSolver.jl's plain rectangular-block collapse —
not a true geometry match). Landed on branch `feat-granular-collapse-example` (off
`mimic-implementation-logic-of-ample`'s tip):

- **Renamed the existing `collapse_problem`** (actually a 1-D elastic self-weight
  column convergence test, not a granular collapse — no plasticity ever exercised)
  out of the way: `collapse.jl`→`column.jl`, `collapse_problem`→`column_problem`,
  `get_collapse.jl`→`get_column.jl`, `get_collapse`→`get_column`,
  `test_collapse.jl`→`test_column.jl`. Verified: 4/4 convergence cases still pass
  unchanged post-rename.
- **New `collapse_problem(nel; kwargs...)`** (`src/home/script/example/collapse.jl`) —
  a real dry granular Drucker-Prager block-collapse example, modeled directly on
  MaterialPointSolver.jl's `2d_collapse.jl` (domain, block size, and material constants
  — `ρ0=2650`, `E` from `Ks=7e5`, `ν=0.3`, `ϕ=19.8°`, cohesionless `c0=0` — all
  hardcoded to match it). Per user's explicit request, `nel` is the *only* argument
  meant to normally change; everything else is a keyword override. Geometry
  (`get_collapse`, `src/home/init/mpts/get_collapse.jl`) reuses the shared
  `mpts_populate` candidate-grid generator and just masks it to a `[0,w0]×[0,h0]` box —
  no new grid-generation math needed.
- **Real bug found and fixed during smoke-testing, user-diagnosed**: the package
  default boundary condition (`:roller` on every side, including the base — normal-
  component-only slip) let particles slide/eject unrealistically at the base. User's
  diagnosis: "problem might be boundary condition. Normally, it is sticky at the
  bottom." Fixed by defaulting `collapse_problem`'s base boundary to `:fixed`
  (sticky/no-slip), matching MaterialPointSolver.jl's own BC treatment (fixes both
  velocity components at the base, only the normal component at the side walls).
  Visually confirmed via before/after scatter plots: fixed the worst artifacts (one
  particle ejected to y=0.14m, others past x=0.7m); some minor stray-particle noise
  remained at coarse resolution with `locking=false` (used only for smoke-test speed,
  not the intended default) — attributed to the same grid-crossing instability
  CLAUDE.md's Known bugs section already documents, not a new bug.
- Added a `"collapsing"` showcase branch to `welcome_log` (`src/home/api/io/logs.jl`),
  mirroring the existing `"slumping"`/`"heating"` pattern.
- Verified: clean load, `test_column.jl` 4/4 (rename didn't break it),
  `test_workflow.jl`'s 96-case sweep unchanged at the documented 71 passed/25 failed
  baseline (this change doesn't touch anything that sweep exercises). CLAUDE.md updated
  throughout (Known bugs' `test_collapse.jl`/`test_column.jl` entry, "Operating the
  package"'s example list with a new `collapse_problem` bullet).
- **Not done in this pass**: a same-geometry timing/visual comparison against
  MaterialPointSolver.jl using the new `collapse_problem` (the original motivation) —
  user said "can i try?" and took over running it themselves partway through; no
  re-comparison numbers were generated by Claude after the BC fix landed. Natural next
  step if picked up again.
- **Also reverted before commit**: the user had experimentally flipped
  `get_default().nonloc.status` from `true`→`false` in `defaults.jl` directly (own
  edit, not Claude's) while investigating non-local regularization's cost earlier in
  the session — explicitly confirmed as "just an experiment... revert when committing",
  so `defaults.jl` was restored to its original `status=true` before this commit.

Scripts used (scratch, not committed): `mps_bench/smoke_test.jl` and
`mps_bench/matched_run.jl` in the session scratchpad, `ep_bench.jl` alongside them.
Used `basis.which="bsmpm"` on the ElastoPlasm side specifically to dodge the known
`gimpm`/grid-crossing `DomainError` (CLAUDE.md Known bugs) — hit it once with
`gimpm` at this particle count before switching.

**Why the archive repo, not current `main`**: investigation found current
`main`/v0.5.1 of `MaterialPointSolver.jl` was refactored into a minimal low-level
core with no runnable example anywhere in its `test/`/docs for that API version
(confirmed via a dedicated Explore agent, not assumed) — resolved by using
`LandslideSIM/Archive_MaterialPointSolver.jl_paper` instead, per the user's own
redirect ("maybe here" → the LandslideSIM org) after being asked to choose between
that, writing a from-scratch kernel loop, or dropping the comparison. `:CPU` was
independently confirmed to be a real, fully-implemented first-class device option in
that vendored solver source (not a fallback) before committing to the plan.

## Queued, not started (superseded — was)

**Performance comparison vs. `MaterialPointSolver.jl`** (github.com/ZenanH/
MaterialPointSolver.jl / LandslideSIM org) — user's idea, explicitly deferred:
"We will tackle this later this day... I will instruct you to do comparison later on."
Also Julia + KernelAbstractions.jl (same backend-dispatch mechanism as ElastoPlasm);
Emmanuel Wyser (this repo's git user) co-authored their Computers & Geotechnics paper,
and their README lists elastoPlasm.jl as a sibling package — not an arbitrary
third-party comparison. Already shallow-cloned to
`c:\Users\lili8\Documents\GitHub\MaterialPointSolver.jl` for inspection.

Closest matched scenario found: their `test/Finite_strain` (3D dry granular collapse,
Drucker-Prager, Hencky finite-strain elasticity, MUSL-style double-mapping P2G/G2P) vs.
ElastoPlasm's `slump_problem`/`collapse_problem`. Blockers to actually running it:
(1) their test file imports `CUDA`/`WGLMakie` and calls `meshbuilder`/`animation` from
two sibling packages (`MaterialPointGenerator`, `MaterialPointVisualizer`) not in their
`Project.toml` — need adding separately; (2) `init_dev` needs switching from `:cuda` to
`:cpu` for a fair CPU-vs-CPU comparison (untested whether their `:cpu` backend path is
actually exercised/working, since KernelAbstractions ports sometimes only get tested on
GPU in practice).

Planned approach when resumed: write a standalone 2D script against their public API
(`init_conf`/`init_grid`/`init_mpts`/`mpmsolver!`) mirroring `slump_problem`'s geometry
(`L=[64.1584, 64.1584/4.0]`, `nel=[40,10]`) and material constants (`E`/`ν`/`ρ`/`ϕ`),
`strain=finite`, `dev=:cpu`, then time it against `elastoplasm!(jld2;
workflows=[elastodynamic!])` for matched total simulated time. Rough estimate given to
the user: ~1–1.5h wall-clock, mostly unattended precompile (`CUDA.jl`/`WGLMakie` are
heavy even in CPU mode) — biggest risk is their `:cpu` backend being untested/broken,
which would turn into real debugging time.

## Standing session notes

- User prefers `function ... end` block syntax over compact `f(...) = expr` Julia
  syntax (saved as its own memory: `julia-function-syntax-style`).
- User prefers saving ongoing-task-state traces in this repo's own `.claude/` folder
  in addition to the cross-session global memory (saved as its own memory, on the
  `track-ongoing-task-state` file).
- Watch for subagents spawning further subagents instead of doing delegated work
  directly and reporting real numbers — happened repeatedly earlier this session;
  redirect explicitly and check `ListAgents` before trusting a "done" report.
