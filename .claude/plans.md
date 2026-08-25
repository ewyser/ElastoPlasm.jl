# Queued plans

Local, gitignored, durable across sessions — unlike the global plan-mode scratch file
(`.claude/plans/*.md` under the user's home directory), which gets reused/overwritten
by unrelated future `/plan` invocations. Check this file at the start of a session on
this repo, alongside `.claude/task-state.md`.

---

# Dispatch transfer scheme via a `Basis` type parameter

**Status: queued, not started.** Approved design, not yet implemented. Start a new
branch off the current `dev` tip (which includes the merged typed-tensor work) when
picking this up — e.g. `feat-dispatch-transfer-scheme-on-basis`.

## Context

Currently `std_p2n`/`tpic_p2n`/`apic_p2n` are three separately-named `@kernel` functions,
selected once at solver-construction time by `init_mapsto`'s `if/elseif` on
`instr.transfer.trsfr` into the `Cairn` closure table (`src/home/core/solver/explicit/
mapsto/mapsto.jl:13-36`), then invoked uniformly as `solver.cairn.mapsto.map.p2n!(...)`.
The goal: replace this with ordinary Julia multiple dispatch — give `Basis` a
`transfer::TR<:AbstractTransfer` field/type parameter (mirroring its existing
`kind::K<:AbstractBasis`, which already dispatches `eval_basis`/`shpfun!` onto
`BSplineBasis`/`GimpBasis`/`LinearBasis`/`MLSBasis`), and rewrite the three `p2n!`
kernels as one function name with three methods dispatching on `Basis{...,TR}`.

**Confirmed feasible and consistent with an existing precedent in this exact file**:
`shpfun!` (`basis/mlsmpm.jl`) is already a `@kernel`-decorated function with multiple
methods dispatching on `Basis`'s `K` type parameter, proving `@kernel` + type dispatch
coexist cleanly here — the same mechanism this plan extends to `TR`.

**Soundness**: CLAUDE.md (line 80-83) states `basis.which`/`transfer.trsfr` are "fully
orthogonal — any basis kind pairs with any transfer scheme." That orthogonality is
about combinability (the full `(K,TR)` cross-product stays constructible either way),
not about which struct's type parameters carry the tag — placing `TR` on `Basis`
doesn't collapse that cross-product. `Basis` already owns "how a particle talks to the
grid" (topology, shape-function kind, cached `N`/`∂N`); transfer scheme is the other
half of the same question — *what* gets transferred once you know *which nodes*.
Putting both on `Basis` is a reasonable design call, not just a tolerated one.

**Confirmed, not assumed**: thermal and mechanical `mapsto` share one `init_mapsto`
call (`get_solver.jl:35`) and thus one `Cairn` `:p2n!` slot — so setting
`transfer.trsfr="apic"`/`"tpic"` with a thermal-only workflow **already throws
`MethodError` today** (no thermal-phase method exists on `apic_p2n`/`tpic_p2n`). This
plan preserves that exact behavior — a thermal call with `TR<:ApicTransfer`/
`TpicTransfer` will still throw `MethodError` (no thermal method for those `TR`s),
just dispatched via `Basis`'s type instead of via Cairn-then-argument-type. No new
design needed there; not an improvement or a regression, a like-for-like preservation.

## Implementation

**1. New marker types** — in `src/boot/needs/types/abstract.jl`, alongside
`AbstractBasis`: add `AbstractTransfer`. In `src/boot/needs/types/concrete/basis/
basis.jl` (or a new sibling file if it reads cleaner — match whichever convention the
`AbstractBasis` concrete kinds use, check `bsmpm.jl`/`gimpm.jl` for the exact pattern
before choosing): add `StdTransfer`/`TpicTransfer`/`ApicTransfer <: AbstractTransfer`,
plain marker structs (no fields needed — the transfer schemes carry no per-scheme
config data unlike `GimpBasis`'s `stencils`).

**2. `get_transfer(which::String) -> AbstractTransfer`** — mirrors `get_basis`
(`basis.jl:23-35`) exactly: `"std"→StdTransfer()`, `"tpic"→TpicTransfer()`,
`"apic"→ApicTransfer()`, else `error`.

**3. `Basis` struct** (`basis.jl:56-67`): add `TR<:AbstractTransfer` as a **trailing**
type parameter (`Basis{T1,T2,D,NN,K<:AbstractBasis,TR<:AbstractTransfer}`) and a
`transfer::TR` field, placed after `kind::K`. Trailing, so the ~13 existing `K`-only
dispatch sites (`Basis{T1,T2,D,NN,K}` pattern in `eval_basis`/`shpfun!`/etc.) need no
edit — Julia leaves an unwritten trailing parameter free, the same mechanism already
verified this session for `Point`'s `ST`/`SC`/`SK`.

**4. `setup_basis.jl`**: call `get_transfer(solver.transfer.trsfr)`, pass the instance
as the new positional field argument, and add `typeof(transfer)` to the `Basis{...}`
constructor's type-parameter list (`setup_basis.jl:40`).

**5. Rewrite `std_p2n`/`tpic_p2n`/`apic_p2n` → one `p2n!`** across `mapsto/fun/p2n/
{std,apic,tpic}.jl`: same function name `p2n!` in all three files, each method's
`basis` argument typed `Basis{T1,T2,D,NN,K,TR}` with `TR<:StdTransfer`/`TpicTransfer`/
`ApicTransfer` added to its `where` clause (`K` stays free/untyped in these, matching
how they don't currently care about basis kind). `std.jl`'s thermal-phase overloads
(different `mesh::MeshThermalPhase` argument, no basis-kind-or-transfer dispatch today)
keep their existing name unless full consistency is wanted — renaming them to `p2n!`
too is purely cosmetic (their dispatch is already on `mesh`'s type, not `basis`'s, so
no `TR` constraint is needed on them), not a functional change.

**6. `Bᵢⱼ!` (currently only registered/called for APIC)**: add trivial no-op methods
for `TR<:StdTransfer`/`TpicTransfer` alongside the real `ApicTransfer` implementation,
so it can be registered and called unconditionally instead of behind the
`if solver.transfer.trsfr == "apic"` check.

**7. `init_mapsto`** (`mapsto.jl:13-36`): replace the `if/elseif` on `instr.transfer.trsfr`
for the `:p2n!`/`:Bᵢⱼ!` Cairn entries with unconditional registration
(`mapsto[:map][:p2n!] = p2n!(CPU())`, `mapsto[:map][:Bᵢⱼ!] = Bᵢⱼ!(CPU())`) — dispatch now
happens via `basis`'s type at call time, so `init_mapsto` no longer needs to know
`instr.transfer.trsfr` at all for these two entries. Keep the `ArgumentError` on an
unrecognized `trsfr` string — that validation still needs to happen somewhere; either
leave a thin check in `init_mapsto` or move it into `get_transfer`'s `else` branch
(already there, per step 2) and drop the now-redundant one from `init_mapsto`.

**8. `mapsto.jl:80`**: replace `if solver.transfer.trsfr == "apic"` with an
unconditional `solver.cairn.mapsto.map.Bᵢⱼ!(mpts,mesh,basis;...)` call, relying on the
new no-op methods from step 6 to make it a correctly-typed no-op for `std`/`tpic`.

## Addendum: `Point`'s `Δnp`/`Bᵢⱼ`/`Dᵢⱼ` fields

Separately flagged: `Point`'s `Δnp`/`Bᵢⱼ`/`Dᵢⱼ` fields (`lagrangian.jl:72-76`, tagged
`# basis-related quantities`/`# APIC-related` in the comments already) belong on
`Basis` instead, consistent with this refactor's direction. Investigated each (grep,
not assumed):

- **`Δnp` is dead code — confirmed, not just suspected.** Zero reads anywhere in
  `src/`; its only occurrences are its own field declaration, its zero-init in
  `setup_mpts.jl:146`, and one docstring mention in `basis.jl` citing it as a
  *precedent example* for the `N`/`∂N` caching pattern (ironic, since the precedent
  itself was never wired up). **Delete it entirely** — same treatment as the `P`/
  `elast::E` dead fields removed earlier this session. Update `basis.jl`'s docstring
  to stop citing it as a precedent (rephrase generically, or cite `N`/`∂N` themselves
  as their own precedent instead).
- **`Bᵢⱼ`/`Dᵢⱼ` are genuinely live and 100% APIC-exclusive already, just not in
  storage.** Confirmed: `Dij_nd` (`shpfun.jl:15-30`, writes `mpts.Dᵢⱼ`) is only
  registered (`ignite.jl:20`) and only called (`ignite.jl:38-40`) when
  `transfer.trsfr == "apic"`; `Bij.jl`'s `Bij` kernel (writes `mpts.Bᵢⱼ`) and
  `apic.jl`'s `p2n!` (reads both) are APIC-only files. So today every `std`/`tpic` run
  still allocates two unused `Vector{SMatrix}` arrays on every `Point` for nothing.
  Move both onto `ApicTransfer` (not `Basis` directly) — `ApicTransfer{T,D}` gains
  `Bᵢⱼ::Vector{SMatrix{D,D,T,D*D}}`/`Dᵢⱼ::Vector{...}` fields, exactly mirroring
  `GimpBasis{T,D,NN}`'s own `stencils::NTuple{D,UnitRange{T}}` field (`gimpm.jl:7-15`)
  — kind-specific data living on the kind struct. Unlike `GimpBasis{T,D}()`'s zero-arg
  constructor, `ApicTransfer` needs `nmp` to size its vectors, which isn't available
  inside `get_transfer` the way `get_basis` works — construct it directly in
  `setup_basis.jl` instead (where `nmp = problem.mpts.nmp` is already in scope), e.g.
  `get_transfer(solver.transfer.trsfr, T2, D, nmp)`, with `StdTransfer()`/
  `TpicTransfer()` staying zero-arg for the other two schemes.
- **Consumers to update**: `shpfun.jl`'s `Dij_nd` and `Bij.jl`'s `Bij` write to
  `basis.transfer.Dᵢⱼ[p]`/`.Bᵢⱼ[p]` instead of `mpts.Dᵢⱼ[p]`/`.Bᵢⱼ[p]`; `apic.jl`'s
  `p2n!` reads the same way. All three already take `basis` as an argument, so this is
  a pure field-access rename, no signature changes needed.
- **`Point`/`PointSolidPhase`/`setup_mpts.jl`**: remove `Δnp`/`Bᵢⱼ`/`Dᵢⱼ` from
  `Point`'s struct body (`lagrangian.jl:66-96`) and their construction
  (`setup_mpts.jl:146,148-149`). `TM` stays in `Point`'s type parameter list
  regardless (still needed to forward into the embedded `PointSolidPhase{...,TM,...}`
  field), just no longer used by a field directly on `Point` itself.
- **Persistence note, real and worth stating plainly**: `Bᵢⱼ`/`Dᵢⱼ` currently persist
  as part of `ic["problem"]` (via `Point`); after this move they persist as part of
  `ic["basis"]` (via `Basis`/`ApicTransfer`) instead. This is a deliberate JLD2 layout
  change, not a side effect to gloss over — call it out in the CLAUDE.md Persistence
  section update, and explicitly round-trip-test it (save/load a `Point`+`Basis` pair
  under `transfer.trsfr="apic"`, confirm `Bᵢⱼ`/`Dᵢⱼ` values and types survive).
- **`ignite.jl`'s `if instr[:transfer][:trsfr] == "apic"` / `if solver.transfer.trsfr ==
  "apic"` guards** (`init_ignite`/`ignite`, gating `Dᵢⱼ!`'s registration/call): can stay
  as string checks (`Dᵢⱼ!` is orthogonal to the `p2n!`/`Bᵢⱼ!` dispatch this plan's main
  body converts) — converting these too is optional polish, not required for the field
  relocation to work, since `Dij_nd`'s only real dependency on knowing "am I APIC" is
  which array it's allowed to write into, which the relocation itself already enforces
  (the field simply doesn't exist under `std`/`tpic`).

## What stays unchanged

- `solver.transfer.trsfr` (the string) stays in the config — still needed by
  `get_transfer` at `setup_basis` time, and by anything reading it for display/logging.
  This plan doesn't remove the config string, only the runtime dispatch mechanism.
- `solver.transfer.musl` (MUSL reprojection) is an orthogonal boolean, untouched —
  stays a runtime `if` check (`mapsto.jl:70`, `init_mapsto`'s musl branch), not folded
  into `TR`.
- `n2p!`/`solve!`/`σᵢ!` and the `augm.*` Cairn entries are untouched — this plan only
  targets the `:p2n!`/`:Bᵢⱼ!` slots, which are the ones actually keyed on transfer
  scheme.

## Verification

1. `using ElastoPlasm` — clean load, no `@warn "Failed to include..."`.
2. Numerical equivalence: for representative particle/mesh state, confirm the new
   dispatched `p2n!` produces bit-identical `mesh.s.m`/`.mv`/`.oobf` output to the old
   `std_p2n`/`tpic_p2n`/`apic_p2n` for each scheme — compare against the pre-change
   code, not just "looks equivalent".
3. Full smoke matrix: `slump_problem` with `transfer.trsfr ∈ {"std","tpic","apic"}` ×
   `plast.constitutive ∈ {"DP","J2"}`, all `success=true`.
4. Confirm the thermal + non-`"std"`-transfer combination still throws the same
   `MethodError` shape as before (per the "what stays unchanged" note above) — don't
   let this regress silently into success or into a different, confusing error.
5. `test_workflow.jl`'s full 96-case sweep — varies `transfer.trsfr` already, so this
   is the most direct regression check available: confirm unchanged at 25/96 failures,
   same per-basis-kind distribution documented in CLAUDE.md.
6. `test_collapse.jl`, `test_basis.jl`.
7. `test_performance.jl` — report before/after; this is the one place a real
   performance difference (closure-table vs. compiled dispatch) might actually show
   up, worth measuring rather than assuming either way.
8. Grep repo-wide for `std_p2n`/`tpic_p2n`/`apic_p2n` afterward to confirm no leftover
   references to the old names.
9. Grep repo-wide for `Δnp` and `mpts.Bᵢⱼ`/`mpts.Dᵢⱼ` afterward — should return zero
   hits (the first deleted, the latter two relocated to `basis.transfer.*`).
10. JLD2 round-trip: save/load a `Point`+`Basis` pair under `transfer.trsfr="apic"`,
    confirm `basis.transfer.Bᵢⱼ`/`.Dᵢⱼ` values and types survive at their new home.
11. Allocation check: confirm a `std`/`tpic` run's `Basis`/`Point` construction no
    longer allocates the (now `ApicTransfer`-only) `Bᵢⱼ`/`Dᵢⱼ` vectors — the concrete
    payoff of moving them off `Point` rather than just leaving them there.

Update CLAUDE.md: the `basis.which`/`transfer.trsfr` orthogonality note (line 80-83)
needs a follow-up explaining `Basis` now also carries `TR`, and why — write it the
same "done, with tradeoffs stated" style CLAUDE.md already uses for e.g. the
`AbstractElasticity` removal.
