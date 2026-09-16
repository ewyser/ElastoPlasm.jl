# `test_workflow.jl`'s full sweep has genuine, reproducible numerical-instability failures confined to `smpm`/`gimpm`

**Status: open** (root cause narrowed but not fully confirmed; `bsmpm`/`mlsmpm` pass cleanly).

Originally a 96-case, 2D-only sweep; the 3D geometry case is now uncommented
too (see "3D conformity check" below), making it 192 cases (96×2D + 96×3D).
2D baseline (re-measured multiple times across major refactors, unchanged
every time): **smpm 12/24 fail, gimpm 13/24 fail, bsmpm 0/24, mlsmpm 0/24 —
25 failures total.** 3D: smpm 12/24, gimpm 6/24, bsmpm 0/24, mlsmpm 0/24 — 18
failures total, same failure signature, not a new instability. Every
failure is a real thrown exception (`DomainError` or `BoundsError`), caught
and logged by the sweep's own `try/@warn` — not a silent `success=false`
path (`elastoplasm`/`elastoplasm!` have no such path). Full failing case
list, for future diffs:
`smpm_{finite,infinitesimal}_{std,tpic,apic}_lock{true,false}_muslfalse` (12)
and
`gimpm_{finite,infinitesimal}_{std,tpic,apic}_locktrue_musl{true,false}`
minus `gimpm_finite_tpic_locktrue_musltrue`, plus
`gimpm_{finite,infinitesimal}_tpic_lockfalse_muslfalse` (13).

## Root-cause theory, narrowed but not fully confirmed

The classic MPM cell/grid-crossing instability. `smpm`'s piecewise-linear
tent shape function has a discontinuous `∂N` at `δx=0` (jumps from `1/h` to
`-1/h` exactly at a mesh node) — a particle crossing an element boundary
sees a discontinuous jump in shape-function gradient and thus in computed
strain-rate/stress at that instant. `bsmpm`'s cubic B-spline is
`C¹`-continuous by construction, and `mlsmpm`'s per-particle moment-matrix
reconstruction likewise avoids the jump — both pass 24/24. `gimpm`'s
particle-domain averaging smooths the discontinuity partially, consistent
with its rarer failure rate. Downstream: once a crossing event jolts a
particle's kinematics, `F`/`detF` can drift enough that `volumetric.jl`'s
`ΔJp` (raises the F-bar-averaged Jacobian ratio to a fractional power with
no sign guard) receives a negative value and throws `DomainError` under
`locking=true`; under `locking=false` the drift can instead push the
particle out of the mesh, throwing `BoundsError` in `p2e`/`p2n`. MUSL
reprojection is basis-kind-agnostic and acts as a general stabilizer that
happens to rescue otherwise-marginal `smpm` runs, consistent with its
dominance in the failure split (`musl=false` accounts for essentially all
`smpm` failures). Ruled out: an initial theory blaming `smpm`'s lack of
`basis.type` boundary correction was wrong — that correction only matters
for *wide*-stencil kernels (`bsmpm`'s cubic B-spline reaches beyond a
particle's own element); `smpm`'s `0:1` stencil is compactly supported
within the mesh by construction, needs no correction. Also ruled out:
cross-case state leakage from running the sweep in one process — the
failures reproduce identically re-run in isolation.

Most promising next diagnostic step: log `mpts.s.ΔFᵢⱼ`/`detF` for a single
particle across the timestep where it crosses element boundaries in a
failing `smpm` case, and confirm a step-change coincides with the crossing.
Real fix is likely a higher-order/damped shape function near crossing events
(out of scope to design here), or at minimum a defensive sign/clamp guard on
`ΔJ^dim` in `volumetric.jl` (doesn't fix the root cause, but closes the
crash).

## Sweep restructuring and nonlocal coverage

`test_workflow.jl`'s config-case generators were restructured to match
actual solver config sections, and nonlocal regularization is now covered by
the sweep. The old `generate_fwrk_cases()` bundled
`deform`/`trsfr`/`locking`/`musl` into one flat, ad-hoc `fwrk` NamedTuple
with no counterpart in the real solver config (`ExplicitSolver` keeps
`strain`/`stab` separate, and `trsfr`/`C_pf` live inside `basis` — no `fwrk`
field anywhere in that struct; a `fwrk` field does exist, but only on the
separate, unused `DynamicRelaxationSolver`, where it plays a different
role). Replaced with four generators, one per real config section —
`generate_strain_cases()`, `generate_transfer_cases()`,
`generate_stab_cases()`, `generate_nonloc_cases()` — each returning
NamedTuples shaped exactly as `slump_problem` expects, mirroring how
`generate_basis_cases()` already worked. The main loop crosses these four
directly and passes them straight through, with no manual field-unpacking.
`generate_nonloc_cases()` is new: nonlocal regularization (`nonloc.status`)
had no sweep coverage at all before this — every prior run of
`test_workflow.jl` hardcoded `nonloc.status=false`, despite nonlocal having
had a real, long-standing correctness bug (see
`bug/fixed/nonlocal-regularization-on2-and-asymmetry-bug.md`).

Re-ran the (now 2D-only again — 3D re-commented out locally, see "3D
conformity check" below) sweep with `nonloc` as a genuine axis: **141/192
passed.** `smpm`'s 12 failing configs are unchanged by nonlocal (just doubled
to 24/48 by the new axis — nonlocal doesn't interact with that instability
either way); `bsmpm`/`mlsmpm` stay fully clean (48/48) regardless of
`nonloc.status`. `gimpm` goes from 13/24 (baseline) to 21/48 — the expected
doubling (26) **plus exactly one new failure**: `gimpm, finite, tpic,
locking=true, musl=true` passes with `nonloc.status=false` but fails with
`nonloc.status=true` (`BoundsError: attempt to access 400-element
Vector{...} at index [16346]` in `element_to_nodes_topology`, `tplgy.jl:15`
— a particle drifting far outside the mesh, same failure mode as the
grid-crossing instability above, just newly triggered on an otherwise-stable
configuration). This narrow, `gimpm`+`tpic`-specific case is not
root-caused — see
`bug/known/test-workflow-nonlocal-gimpm-tpic-boundserror.md`, since
nonlocal regularization's plastic-strain averaging could plausibly be
feeding back into the same kind of drift.

## 3D conformity check — done, via `test_workflow.jl` itself

The 3D geometry case in `generate_geometry_cases()`
(`L=[64.1584,64.1584/4,64.1584/4]`, `nel=[40,10,10]` — matching the 2D case's
resolution exactly, not a coarser stand-in) is now **uncommented and part of
the regular 192-case sweep** (96×2D + 96×3D). Full run: **149/192 passed
(71/96 2D + 78/96 3D)**. 2D reproduces the documented baseline exactly (smpm
12/24, gimpm 13/24, bsmpm 0/24, mlsmpm 0/24 fail). 3D: **smpm 12/24, gimpm
6/24, bsmpm 0/24, mlsmpm 0/24 fail** — same shape as 2D, and *better* than 2D
on `gimpm` specifically (6 fail vs. 13) rather than worse. Every 3D failure
fits the already-documented 2D root cause: mostly the same `DomainError` (a
negative value into `volumetric.jl`'s unguarded fractional-power `ΔJ^dim`,
e.g. `-123.7`, `-567.6`), plus 3 `BoundsError`s (`attempt to access
4000-element Vector{...} at index [-303]` etc.) matching the documented
"particle drifts out of mesh under `locking=false`" failure mode.
**Conclusion: 3D is not meaningfully less stable than 2D** for the
basis-kind/transfer-scheme/locking/musl configuration space this sweep
covers — no new 3D-specific failure category turned up. Plasticity *does*
actually run in this sweep in both dimensions, despite `plast.status` staying
at its `false` default — see `bug/known/plast-status-dead-config-flag.md`:
`elastoplast()` (invoked by the `elastoplastic!` workflow, which every case
here runs) calls `retmap!` unconditionally, with no `if solver.plast.status`
gate at that or any call site in `src/home/` — `plast.status` is read nowhere
except a docstring comment. Which workflow function you pass to
`elastoplasm!`/`elastoplasm` (`elasto` via `elastodynamic!`, vs.
`elastoplast` via `elastoplastic!`) is what actually decides whether the
plastic corrector runs, not the config flag.

**Could not reproduce the previously-reported "3D APIC hits `DomainError`
even with `stab.locking=false`" issue** (see the thermal-solution section of
`.claude/docs/planned-improvements.md`) under the current code, across
several attempts: plain `slump_problem` 3D+APIC (`bsmpm`, `locking=false`),
same with `plast.status=true` enabled, and `thermal_problem` 3D+APIC — all
ran clean (`success=true`). Either it was fixed as a side effect of unrelated
later work, or it needs a more specific trigger (finer resolution,
`collapse_problem`-style violent impact, a different `nel`/domain) not yet
found. Leaving this as unconfirmed/unreproduced rather than deleting it,
since "couldn't reproduce" isn't the same as "fixed."

Found instead, a real and separately-fixed bug along the way: see
`bug/fixed/thermal-problem-2d-hardcoded-bcs.md`.
