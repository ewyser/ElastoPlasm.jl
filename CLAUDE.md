# ElastoPlasm.jl — Notes for Claude

Working notes on this repo's architecture and conventions, kept up to date as the
codebase evolves. Read this first when picking up work here.

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
implements its concrete subtype(s) (no central `abstract.jl`) — see "Core types" below
for the directory layout. Adding a new file under `types/` means adding it to that
explicit list in `boot.jl`, in the right dependency position; it is NOT picked up
automatically the way `home/` files are.

## Core types

- `Mesh{T1,T2,D}` — Eulerian background grid. Carries no connectivity (`e2n`/`e2e` live
  on `Basis`).
- `Point{T1,T2,D,CM<:AbstractConstitutiveModel,TM,TV,TS,ST<:AbstractStrain,SC<:AbstractStress,SK<:AbstractStress}`
  — material points (Lagrangian). `ST`/`SC`/`SK` are the typed strain/stress storage on
  `mpts.s` (see "Typed strain/stress tensor storage" under Planned improvements);
  `mpts.s.σᵢⱼ[p]` returns a `CauchyStress`, **not** an `SVector` — read it with
  `get_voigt(...)`. Carries no `NN` or connectivity (those live on `Basis`).
  `mpts.x :: Vector{SVector{D,T2}}` — NOT a matrix; index with `getindex.(x, i)` or
  iterate, never `x[i, :]`. `mpts.s.cmp::Vector{CM}` is the per-particle
  constitutive-model bundle (see "Typed constitutive-model abstraction"). `CM` resolves
  to `DruckerPrager` or `VonMises` depending on `plast.constitutive` (see "DP/J2 retmap
  kernel unification"); `ST` to `LogarithmicStrain`/`InfinitesimalStrain` depending on
  `strain.deform`.
- `MechanicalProblem{T1,T2,D,CM,TM,TV,TS,ST,SC,SK} <: AbstractProblem{T1,T2,D,SP}`
  (`src/boot/needs/types/problem/problem.jl`) — bundles `mesh::Mesh`+
  `mpts::Point`+`time::Time` as the IC-defining part of a simulation, built via
  `setup_problem` (see "`Problem` type decoupling..." under Planned improvements).
  Deliberately excludes `Basis`/`Solver` — those stay independent types.
- `Basis{T1,T2,D,NN,K<:AbstractBasis,TR<:AbstractTransfer}`
  (`src/boot/needs/types/basis/basis.jl`) — independent struct owning ALL
  mesh/point connectivity (`e2n`, `e2e`, `p2n`, `p2e` — `e2e` also doubles as the
  neighbor-search structure for nonlocal plastic-strain regularization, see Known bugs),
  plus `kind::K`
  (dispatches shape-function evaluation, `eval_basis`/`shpfun!`), `transfer::TR`
  (dispatches P2G/G2P transfer-scheme kernels, `p2n!`/`Bij` — see "Transfer scheme
  dispatch" below), `NN` (nodes per element — only known once mesh+basis kind are fixed,
  hence its own type param), `type` (per-axis node boundary-layer classification,
  `BSplineBasis`-specific), and the shape-function cache `N`/`∂N` (see below).
- Construction order is always **mesh → material points → basis** (basis depends on
  both). Function argument order is always **`(mpts, mesh, basis, ...)`**.
- Shape-function evaluation is `eval_basis(mpts, mesh, basis, ip, nn)`. Concrete kinds
  live in `basis/{bsmpm,gimpm,smpm,mlsmpm}.jl`, each defining its own
  `stencils::NTuple{D,UnitRange{T}}` field on the kind struct.
- **Shape values/gradients are cached, not recomputed per call site.** `Basis.N` /
  `Basis.∂N` (plain `Matrix{T2}` / `Array{T2,3}`, sized `(NN,nmp)` / `(NN,D,nmp)` — see
  "Precision/StaticArrays gotchas" below for why not `Vector{SVector}`/`Vector{SMatrix}`)
  are populated once per particle per timestep by the `shpfun!` ignite kernel
  (`src/home/core/common/shpfun.jl`). Every P2G/G2P/update kernel reads `basis.N[nn,p]` /
  `∂Nrow(basis.∂N,nn,p,Val(D))` (a small helper in `basis.jl` extracting a gradient row
  via a `Val`-unrolled `ntuple`, since `M[nn,:]` on an `SMatrix` allocates for a runtime
  `nn`) instead of calling `eval_basis` inline. `dynamic_relaxation`'s ~14 call sites are
  the one exception — no ignite phase there, so they call `eval_basis` directly.
- `MLSBasis` (`basis/mlsmpm.jl`) is a fourth kind implementing structured-grid Moving
  Least Squares shape functions (Cao et al. 2025, *Comput. Mech.* 75:655–678, §2.3.1 —
  only the structured-grid case). Needs no boundary node-type correction (`basis.type`)
  — a per-particle moment-matrix projection restores partition-of-unity everywhere on
  its own — but needs a once-per-particle `(D+1)×(D+1)` matrix build+inversion, which is
  why the `shpfun!` caching above exists. `basis.which` and `basis.trsfr` are fully
  orthogonal — any basis kind pairs with any transfer scheme, including `mlsmpm`+`std`.

**Transfer scheme dispatch via `Basis`'s `TR` type parameter.** `p2n!` (P2G projection)
and `Bij` (APIC's affine-velocity update) are dispatched on `Basis`'s `TR`
(`StdTransfer`/`TpicTransfer`/`ApicTransfer`), the exact same "config string → concrete
marker instance → multiple dispatch" pattern `kind`/`get_basis` already uses for basis
kind, and the same pattern the DP/J2 `retmap` unification uses for `Point`'s `CM`+`ST`
(see "DP/J2 retmap kernel unification" below). `get_transfer(basis.trsfr, T2, D, nmp)`
(`transfer.jl`) maps `"std"`/`"tpic"`/`"apic"` to `StdTransfer()`/`TpicTransfer()`/
`ApicTransfer{T2,D,L}(nmp)`; `init_mapsto` registers `p2n!(CPU())`/`Bij(CPU())`
unconditionally, no string branch. `Bij`'s `TR<:StdTransfer`/`TpicTransfer` methods are
no-ops, so `mapsto.jl` calls it unconditionally too. Unrecognized `basis.trsfr` strings
fail fast in `get_transfer` at `setup_basis` time.

`ApicTransfer{T,D,L}` carries `Bᵢⱼ`/`Dᵢⱼ` (`Vector{SMatrix{D,D,T,L}}` — `L` is a separate
type parameter, computed as `D*D` only once `D` is concrete inside the inner
constructor, since `D*D` on a bare `TypeVar` isn't valid in a struct's own field-type
declaration). These fields used to live on `Point` unconditionally (`Δnp`/`Bᵢⱼ`/`Dᵢⱼ`,
100% APIC-exclusive, allocated on every particle even under `std`/`tpic`); `Δnp` was
dead code and got deleted, `Bᵢⱼ`/`Dᵢⱼ` moved onto `ApicTransfer` (`basis.transfer.Bᵢⱼ`,
not `mpts.Bᵢⱼ`). They persist as part of `ic["basis"]` now, not `ic["problem"]`. Removing
this dead/wasted per-`Point` allocation measured **-18.3% memory / -5.4% allocs / 0.0%
time** on `test_performance.jl` vs. the pre-refactor baseline.

**Config note**: `basis.trsfr`/`.C_pf` (transfer scheme + PIC/FLIP blend) live under the
`basis` config section (mirroring `Basis` owning both `kind` and `transfer`); `stab.musl`
(MUSL velocity reprojection) lives under `stab`, since it's a stabilization technique
applied regardless of which transfer scheme ran, not transfer-scheme-specific. There is
no `transfer` config section anymore. See "Controlling solver behaviour" below for the
full config layout and a real shallow-merge bug this reorganization surfaced twice.

**Known, deliberate behavior change**: thermal-only workflows now work under any
`basis.trsfr` (previously threw `MethodError` for `apic`/`tpic`) — that restriction was
an accident of `std_p2n`/`tpic_p2n`/`apic_p2n` being separately-named kernels (thermal
only ever had a method under the `std_p2n` name); now that every `p2n!` method shares
one name, thermal's method (which carries no `TR` constraint, since it never needed one)
matches regardless of `TR`. Confirmed via `git worktree` diff against the pre-refactor
commit to be the *only* behavior change from this refactor; not a regression.

## Explicit vs dynamic_relaxation solvers

- `src/home/core/common/` — kernels/wrappers shared by **both** solver paths:
  `tplgy.jl` (`p2e2n`, topology), `shpfun.jl` (`shpfun!`/`Dij_nd`), `ignite.jl`
  (`init_ignite`/`ignite()`). Both `elastodynamic!`/`elastoplastic!` *and*
  `elastoquasistatic!` call the same `ignite(mpts,mesh,basis,solver::AbstractSolver)`
  wrapper, dispatching on the *abstract* solver type since `dynamic_relaxation` doesn't
  require a particular concrete `solution` value (unlike `elastodynamic!`/
  `elastoplastic!`, typed to `ExplicitSolver` specifically — see `solution` in
  "Controlling solver behaviour" below).
- `src/home/core/solver/explicit/` — the primary, actively-used path
  (`elastodynamic!`, `elastoplastic!`). Kernel dispatch table lives on
  `solver.cairn` (a `NamedTuple`, dot-accessed: `solver.cairn.ignite.tplgy!`,
  `solver.cairn.mapsto.map.p2n!`, `solver.cairn.update.deform!`, ...), built by
  `init_ignite`/`init_mapsto`/`init_update`/`init_implicit` in `get_solver.jl`. Every
  kernel/function in this path takes `basis::Basis` as its 3rd positional arg after
  `(mpts, mesh, ...)`.
- `src/home/core/solver/dynamic_relaxation/` — quasi-static path (`elastoquasistatic!`,
  plus the not-fully-wired u-P variant `elastoquasistaticuP!`/`elastouP!`). Runs
  end-to-end through the same `slump_problem`/`elastoplasm(jld2; workflows=[...])`
  pipeline as the explicit solver — but is much slower to smoke-test (a full
  `elastoquasistatic!` run takes ~1 minute vs. seconds for `elastodynamic!`/
  `elastoplastic!`), so budget for that when verifying a change touches this path.

## Known bugs

Bugs discovered but not yet fixed, kept here so they don't get rediscovered from scratch.
**When picking one up: create a new branch off the current branch, named so the fix is
identifiable (e.g. `fix-volumetric-locking-zero-mass-guard`), rather than fixing in place
on an unrelated branch.**

- ~~`test_workflow.jl`'s 96-case sweep reported 47/49 instead of the documented
  71/25 baseline~~ — **fixed**. Root cause: **every one of `mlsmpm`'s 24 cases was
  throwing `MethodError: ... is ambiguous`** — a genuine `shpfun!` dispatch ambiguity
  between the generic method (`src/home/core/common/shpfun.jl`) and `mlsmpm.jl`'s
  MLS-specific method, introduced by a commit that narrowed the generic method's
  `Point` pattern from `Point{T1,T2,D,CM}` (4 explicit type params) to `Point{T1,T2,D}`
  (3 explicit). Both patterns match the exact same set of concrete `Point` types (`CM`
  unconstrained either way) — but Julia's method-specificity algorithm empirically
  treats an explicitly-named-but-free type parameter differently from an implicit
  "rest" pattern when compared against a sibling method (`mlsmpm.jl`'s) constraining a
  *different* parameter (`Basis{...,K<:MLSBasis}`) — a real, obscure dispatch quirk,
  not reviewable by eye (confirmed via `Test.detect_ambiguities` on an isolated
  reproduction). Fixed originally by naming `CM` explicitly again; later made
  structurally robust by constraining the generic method's `K` to
  `Union{BSplineBasis,GimpBasis,LinearBasis}` (explicitly excluding `MLSBasis`) instead
  of relying on the naming quirk — the two `shpfun!` methods' `Basis` type patterns are
  now disjoint by construction. Verified: sweep exactly matches the documented 71/25
  baseline, identical per-basis-kind breakdown (smpm 12 fail, gimpm 13 fail, bsmpm 0
  fail, mlsmpm 0 fail — see the grid-crossing entry below for that baseline's own
  history); zero `shpfun!` ambiguities via `Test.detect_ambiguities`.

  **Methodology lesson from chasing this**: an earlier pass at this investigation
  wrongly concluded `elastodynamic!`/`elastoplastic!` produce exactly zero
  velocity/stress for any `slump_problem` run, treating it as a separate severe bug.
  That was **not a real bug** — it was caused by calling `elastoplasm(...)` (no `!`)
  and then re-loading `mpts` from the JLD2 file afterward. `elastoplasm` opens the file
  read-only and never writes the post-workflow state back — only `elastoplasm!` does —
  so the re-load always returned the untouched, zero-initialized pre-simulation state.
  **Always use `elastoplasm!`, never `elastoplasm`, when the goal is to inspect
  post-run `mpts`/`mesh` state from the JLD2 file afterward** — `elastoplasm`'s return
  value doesn't carry the mutated objects either (just a bare `(;simulation,success)`),
  so there is no way to recover real results from it short of re-reading the file.

- ~~`basis.which` ∈ {`"gimpm"`, `"mlsmpm"`} combined with `stab.locking=true` crashes
  with `InexactError: trunc(Int64, NaN)`~~ — **fixed**: `volumetric.jl`'s `ΔJp` now has
  the same `iszero(mesh.s.m[no])` guard the velocity kernels already had.
- ~~`mpts.s.Δλ` never written under `plast.constitutive="J2"`, silently disabling
  non-local regularization~~ — **fixed**. `J2.jl`'s return-map kernels computed `Δλ` as
  a purely loop-local quantity and never assigned `mpts.s.Δλ[p]`; since `nonlocal.jl`
  gates all three of its branches on `mpts.s.Δλ[p] != 0`, non-local regularization
  never activated under J2. Fix: accumulate `Δλ` across the CPA iterations and store
  the sum. Measured before→after on a 591-particle run: `Δλ != 0` on 0→424 particles;
  DP runs were never affected.
- ~~`finite_DP` never reset `mpts.s.Δλ[p]` on an elastic step~~ — **fixed**: a particle
  that yielded once stayed permanently flagged "yielding" for `nonlocal.jl`'s gate.
  Same one-line reset its infinitesimal sibling already had.
- ~~`_kirchoff_stress`'s 3D Voigt shear ordering disagreed with every consumer~~ —
  **fixed** as a side effect of the tensor port: the old 3D stress builder used
  slots 4/5/6 = xy/xz/yz, but every `p2n!` kernel and `Del`-facing code uses
  4/5/6 = yz/xz/xy. `get_voigt` now uses the consumer convention throughout. Note this
  means **3D finite-strain results change** (2D unaffected — one shear slot only); not
  verified against a reference solution, only made internally self-consistent. 3D
  remains generally undertested (see "3D conformity check" below).
- **`dynamic_relaxation`'s u-P path throws if exercised.** `implicit.jl`'s
  `pt_solve_uP!` calls `instr.cairn.implicit.fint_p2n!(...)`, a key `init_implicit`
  never registers. `elastoquasistatic!` itself doesn't reach this code path, so it
  isn't blocked by this; out of scope until someone designs what `fint_p2n!` should do.
- **`dynamic_relaxation`'s plain (non-u-P) path only actually supports
  `strain.deform="finite"`**, despite a comment in `implicit.jl` claiming the opposite.
  `pt_solve!`/`relax` always call `oobf_assembly!` with the 4-trailing-arg
  `ST<:LogarithmicStrain` (finite-strain) shape; under `"infinitesimal"` only a 3-arg
  method exists, throwing `MethodError`. Confirmed pre-existing, not caused by any
  recent refactor. `test_column.jl` and every other verified `dynamic_relaxation` run
  use the default `"finite"`, which is why this went unnoticed. Someone needs to either
  fix `oobf_assembly`'s infinitesimal method to actually get called, or make
  `"infinitesimal"` genuinely unsupported by construction (an explicit error) rather
  than an opaque `MethodError` deep in a kernel.
- ~~`test_collapse.jl` fully non-functional~~ — **fixed, and since renamed**.
  `collapse_problem` (renamed from `ic_collapse`) was rewritten to the current
  `get_solver → setup_geometry → setup_problem → setup_basis → export_problem`
  pipeline; `test_collapse.jl` updated to read `file["ic/problem"]`/`file["cfg/solver"]`
  instead of the old `ic["mesh"]`/`ic["mpts"]`/`ic["cmpr"]`/`cfg["instr"]` layout. All 4
  `nel` convergence cases pass, including `@test errors[k+1] < errors[k]` — a real
  physics-correctness signal. Also the first real exercise of `dynamic_relaxation`'s
  `cmp`-rerouted reads (see "Typed constitutive-model abstraction" below), now
  confirmed working. **Since renamed**: this test is a 1-D elastic self-weight column
  convergence check (no plasticity ever exercised), not a granular collapse — the
  misleading name was freed up for an actual granular-collapse example, `collapse_problem`
  (see the `collapse_problem`/`column_problem` bullets under "Operating the package"
  below). `collapse_problem`→`column_problem`,
  `collapse.jl`→`column.jl`, `get_collapse`→`get_column`,
  `test_collapse.jl`→`test_column.jl`; the 4/4 convergence result above is unchanged
  post-rename.
- **`test_workflow.jl`'s full 96-case sweep has genuine, reproducible
  numerical-instability failures confined to `smpm`/`gimpm` — `bsmpm`/`mlsmpm` pass
  cleanly.** Current baseline (re-measured multiple times across major refactors,
  unchanged every time): **smpm 12/24 fail, gimpm 13/24 fail, bsmpm 0/24, mlsmpm
  0/24 — 25 failures total.** Every failure is a real thrown exception
  (`DomainError` or `BoundsError`), caught and logged by the sweep's own `try/@warn` —
  not a silent `success=false` path (`elastoplasm`/`elastoplasm!` have no such path).
  Full failing case list, for future diffs:
  `smpm_{finite,infinitesimal}_{std,tpic,apic}_lock{true,false}_muslfalse` (12) and
  `gimpm_{finite,infinitesimal}_{std,tpic,apic}_locktrue_musl{true,false}` minus
  `gimpm_finite_tpic_locktrue_musltrue`, plus
  `gimpm_{finite,infinitesimal}_tpic_lockfalse_muslfalse` (13).

  **Root-cause theory, narrowed but not fully confirmed**: the classic MPM
  cell/grid-crossing instability. `smpm`'s piecewise-linear tent shape function has a
  discontinuous `∂N` at `δx=0` (jumps from `1/h` to `-1/h` exactly at a mesh node) — a
  particle crossing an element boundary sees a discontinuous jump in shape-function
  gradient and thus in computed strain-rate/stress at that instant. `bsmpm`'s cubic
  B-spline is `C¹`-continuous by construction, and `mlsmpm`'s per-particle
  moment-matrix reconstruction likewise avoids the jump — both pass 24/24. `gimpm`'s
  particle-domain averaging smooths the discontinuity partially, consistent with its
  rarer failure rate. Downstream: once a crossing event jolts a particle's kinematics,
  `F`/`detF` can drift enough that `volumetric.jl`'s `ΔJp` (raises the F-bar-averaged
  Jacobian ratio to a fractional power with no sign guard) receives a negative value
  and throws `DomainError` under `locking=true`; under `locking=false` the drift can
  instead push the particle out of the mesh, throwing `BoundsError` in `p2e`/`p2n`.
  MUSL reprojection is basis-kind-agnostic and acts as a general stabilizer that
  happens to rescue otherwise-marginal `smpm` runs, consistent with its dominance in
  the failure split (`musl=false` accounts for essentially all `smpm` failures).
  Ruled out: an initial theory blaming `smpm`'s lack of `basis.type` boundary
  correction was wrong — that correction only matters for *wide*-stencil kernels
  (`bsmpm`'s cubic B-spline reaches beyond a particle's own element); `smpm`'s `0:1`
  stencil is compactly supported within the mesh by construction, needs no correction.
  Also ruled out: cross-case state leakage from running the sweep in one process — the
  failures reproduce identically re-run in isolation.

  Most promising next diagnostic step: log `mpts.s.ΔFᵢⱼ`/`detF` for a single particle
  across the timestep where it crosses element boundaries in a failing `smpm` case, and
  confirm a step-change coincides with the crossing. Real fix is likely a higher-order/
  damped shape function near crossing events (out of scope to design here), or at
  minimum a defensive sign/clamp guard on `ΔJ^dim` in `volumetric.jl` (doesn't fix the
  root cause, but closes the crash).
- **`ic_collision`'s initial-velocity plot crashes**: `what_plot_field` errors with
  `type NamedTuple has no field data`, because `collision.jl`'s `plot.what` entries are
  hand-built `NamedTuple`s missing the `data`/`label`/`unit`/`scale`/`cb` fields
  `what_plot_field` expects, instead of going through `get_mpts_variable_config()[name]`
  like `slump_problem` does. Pre-existing, not touched by the `ic_collision` rewrite
  that fixed its undefined-`kwargser`/stale-pipeline issues. Workaround: run with
  `plot.status=false` until fixed.
- ~~`MeshSolidPhase.Mᵢⱼ::Matrix{T2}` allocated a full **dense `nno×nno` matrix**
  unconditionally in `setup_mesh.jl`, regardless of solver~~ — **fixed (removed)**.
  O(nodes²) memory: harmless at small mesh sizes (e.g. 6,601 nodes → ~349 MB) but
  `OutOfMemoryError` at moderate resolution (32,841 nodes → ~8.6 GB) — found while
  trying to run `collapse_problem` at a MaterialPointSolver.jl-matched fine resolution
  (`nel=[368,88]`). Confirmed dead (zero reads anywhere in `src/`, only its own
  declaration/construction) before removing — the explicit solver path already uses
  the separate lumped/diagonal `m::Vector{T2}` nodal mass field for everything; `Mᵢⱼ`
  looks like leftover scaffolding for a never-finished consistent-mass-matrix
  formulation. Removed the field from `MeshSolidPhase`
  (`src/boot/needs/types/problem/eulerian.jl`) and its construction in
  `build_solid_mesh_phase` (`src/home/init/setup_mesh.jl`).
- ~~Non-local plastic-strain regularization (`nonlocal.jl`) was `O(nmp²)` in both memory
  and time, every timestep~~ — **fixed, and a real correctness bug found along the way**.
  `Basis.e2p`/`p2p` (`nmp×nel`/`nmp×nmp`) were **unconditionally dense** (`setup_basis.jl`
  built them via `T1.(spzeros(...))`, but the struct fields were plain `Matrix`, so the
  sparse-zeros broadcast got densified on construction) regardless of
  `solver.nonloc.status` — same class of bug as the `Mᵢⱼ` entry above, and just as live:
  `nonloc.status` **defaults to `true`** (`defaults.jl`), independent of `plast.status`,
  so this ran by default in essentially every simulation. Every timestep also allocated a
  fresh `w = spzeros(T2,nmp,nmp)` and wrote into it via scalar `w[p,q]=...` —
  `setindex!` on a `SparseMatrixCSC` at a previously-zero location is `O(nnz)` (shifts
  every higher-index stored nonzero), making this pattern worse than even a dense fill.
  Measured at `nmp=10133`: 2.213 GiB → 986 MiB per `elastodynamic!` run (`plast.status=
  true`, 1137 steps); at `nmp=40413` the fixed version completed in 97.6s/1.433 GiB while
  the pre-fix version was still running after 30+ minutes (not fully timed to completion
  — analytically it's the expected `O(nmp²)` inner-loop-iteration blowup: the old
  `"p->q"`/`"p<-q"` phases scanned **all** `nmp` particles per yielding particle instead
  of using `basis.e2e` — already a genuinely sparse, `ls`-bounded element-neighbor
  structure — to narrow the search).

  Fix: a per-step CSR element→particle bucket list (`build_el2p`, `tplgy.jl`, `O(nmp+nel)`,
  sequential) combined with `basis.e2e` bounds the neighbor search to `O(nmp×k)` (`k` =
  local particle density within `ls`, resolution-independent); `Basis.e2p`/`p2p` and the
  per-step `w` matrix are gone entirely (pairwise weights are cheap enough — one
  `sqrt`+`exp` — to recompute in both passes rather than cache).

  **While verifying this fix bit-for-bit against the old code, found the old code was
  never actually correct**: it deduplicates each unordered pair `(p,q)` via `if
  w[p,q]==0` before recording it into `p2p`, processing particles in ascending index
  order — whichever particle in a pair is visited *first* claims the pair and writes
  `p2p[q,p]`; the *other* particle's later iteration finds `w` already nonzero and never
  writes its own side. Net effect: for a mutual-neighbor pair, only the *lower-indexed*
  particle ends up receiving that neighbor's contribution in the `"p<-q"` averaging
  pass — the *higher-indexed* particle silently loses it. Confirmed directly: with every
  particle yielding, 8086 of 8677 recorded neighbor-pair entries in the old `p2p` were
  asymmetric (present in only one direction) — the regularization has been missing
  roughly half of each particle's true neighbors, systematically, since this code was
  written. The rewrite computes each particle's own neighbor sum independently (no
  cross-particle writes, no dedup guard) and is verified symmetric and correct by
  independent brute-force recomputation of `W`/`ϵpII[2]` for several particles — it does
  **not**, and structurally cannot, reproduce the old asymmetric numbers; that is the
  fix, not a regression.

## Planned improvements (not bugs, just known follow-up work)

- **`Problem` type decoupling mesh/point generation from basis — done, as a thin
  wrapper.** `MechanicalProblem` bundles `Mesh`+`Point`+`Time` — the expensive,
  IC-defining part of a simulation — as one object, reusing the existing
  `MaterialPointPhase` hierarchy rather than a second "phase" concept. `Basis`/`Solver`
  remain separate, un-bundled types, since `basis.which` and `basis.trsfr` vary
  independently. `types/problem/` holds every file that feeds `MechanicalProblem`
  construction — `constitutive.jl`, `geometry.jl`, `eulerian.jl` (Mesh), `tensor.jl`
  (the `AbstractStrain`/`AbstractStress` tensor types `Point` stores), `lagrangian.jl`
  (Point), `time.jl` (Time), `problem.jl` (the bundle itself) — loaded in that explicit
  order by `boot.jl`'s hand-ordered include list (see "Project shape"), not by
  alphabetical file order.

  `setup_problem` (`src/home/init/setup_problem.jl`) is the entry point — calls
  `setup_mesh`/`setup_material_constants`/`setup_mpts`/`setup_time` internally and
  returns a `Problem`; `slump_problem` (renamed from `ic_slump`) uses it.
  `export_problem` writes a single `ic["problem"]` entry instead of separate
  `ic["mesh"]`/`ic["mpts"]`/`ic["time"]` keys — see "Persistence (JLD2)" below.
  `Basis` stays its own top-level `ic["basis"]` key, deliberately independent of
  `Problem`. Kernel/workflow signatures (`workflow!(mpts,mesh,basis,time,solver)`) are
  untouched — `Problem` is, for now, only an outer-API/persistence convenience, not
  threaded through kernels. `setup_basis` takes `(problem::MechanicalProblem, solver)`.
- **Kernel granularity**: `elast!` is the model to follow — decomposed into small
  functions each doing one specific task, rather than one monolithic kernel body.
  `update.jl`'s `elastoplast`/`elasto` dispatch functions (branch on
  `strain.deform`/`basis.how` inline) would benefit from the same split.
- **Typed constitutive-model abstraction — done.** `mpts.s.cmp::Vector{CM} where
  CM<:AbstractConstitutiveModel` (`constitutive.jl`) bundles the *static* elastic+
  plastic material constants (`Gc`, `Kc`, `Del`, `Hp`, `c₀`, `cᵣ`, `ϕ₀`) into one typed
  object per particle, replacing an old global `cmpr::NamedTuple` threaded through ~15
  signatures. Evolving plastic state (`Δλ`, `ϵpII`, `ϵpV`) deliberately stays *outside*
  `cmp`, as flat mutable vectors on `PointSolidPhase`, since it's written every
  timestep and bundling it into an immutable per-particle struct would force a
  full-struct reconstruction on every write.
  - `PerfectlyElastic` is defined but has **no construction path** — dead scaffolding.
    `DruckerPrager`/`VonMises` both **do** have construction paths, branched on
    `plast.constitutive` — see "DP/J2 retmap kernel unification" below.
  - `PointSolidPhase.rheo::R`/`AbstractRheology` — a parallel, never-read earlier
    attempt at the same per-particle constitutive-data problem — was retired entirely,
    since `cmp` already supersedes it.
  - `dynamic_relaxation`'s `cmp`-rerouted reads are confirmed working for the plain
    (non-u-P) path via `test_column.jl`'s convergence sweep, including its
    correctness assertions. The u-P variant remains unverified (and throws immediately
    if exercised regardless — see Known bugs).
- **DP/J2 retmap kernel unification — done.** `retmap/DP.jl`/`retmap/J2.jl`'s four
  separately-named kernels (`finite_DP`/`infinitesimal_DP`/`finite_J2`/
  `infinitesimal_J2`, picked between via a `plast.constitutive`/`strain.deform` string
  branch) are now one kernel name, `retmap`, with four methods dispatched on `Point`'s
  `CM` (`DruckerPrager`/`VonMises`) and `ST` (`LogarithmicStrain`/`InfinitesimalStrain`)
  — mirroring `elast!`'s `ST`-only dispatch, extended to a second axis.
  `init_update` registers `retmap(CPU())` unconditionally, no string branch.

  `VonMises{T2,D,NSTR,L} <: AbstractConstitutiveModel{T2,D}` — same fields as
  `DruckerPrager` minus `ϕ₀` (J2 yield has no pressure dependence). `setup_cmp`
  branches on a `constitutive::String` keyword to build `Vector{DruckerPrager}`
  (`"DP"`) or `Vector{VonMises}` (`"VM"`); an unrecognized string throws immediately at
  setup time. Config string unified to `"VM"` everywhere. Both `DP.jl`/`J2.jl` share a
  per-model return-mapping core (`_druckerprager_return_map`/`_vonmises_return_map`)
  called by both the finite- and infinitesimal-strain `drucker_prager`/`von_mises`
  methods, rather than the infinitesimal kernels hand-duplicating the algebra
  separately — verified the extracted CPA loop matches the original inline version
  bit-for-bit.

  **Gotcha**: `slump_problem` only plots `"phi0"` (friction angle) in its initial-
  condition figure when `plast.constitutive=="DP"` — `VonMises` genuinely has no `ϕ₀`
  field, so plotting it unconditionally used to crash under `"VM"`. Don't fake a value;
  the field just isn't applicable for a pressure-independent yield surface.
- **Typed strain/stress tensor storage — done.** `mpts.s.ϵᵢⱼ`/`.ϵn`/`.σᵢⱼ`/`.σn`/`.τᵢⱼ`
  hold typed tensor objects from `src/boot/needs/types/problem/tensor.jl`, each
  storing a volumetric+deviatoric additive split:
  - `LogarithmicStrain{S,T,L}` / `InfinitesimalStrain{S,T,L}` (`vol::T`,
    `dev::SMatrix{S,S,T,L}`) — `ϵᵢⱼ`/`ϵn`, picked by `solver.strain.deform`.
  - `CauchyStress{S,T,L}` / `KirchhoffStress{S,T,L}` (`p::T` positive in compression,
    `dev::SMatrix{S,S,T,L}`) — `σᵢⱼ`/`σn` are always Cauchy, `τᵢⱼ` always Kirchhoff.

  The abstract supertypes (`AbstractTensor`/`AbstractStrain`/`AbstractStress`) live in
  `types/problem/tensor.jl` alongside their concrete subtypes — `tensor.jl` is
  deliberately included before `types/problem/lagrangian.jl` in `boot.jl`'s explicit
  list (see
  "Project shape"), since `Point`/`PointSolidPhase` there constrain their `ST`/`SC`/`SK`
  type parameters against these abstracts. `PointSolidPhase`/`Point`/`MechanicalProblem`
  carry `ST`/`SC`/`SK` as **trailing** type parameters specifically so pre-existing
  `Point{T1,T2,D}`-pattern kernel signatures needed no edit.

  Two things worth knowing before touching this:
  - **The stored split is representational, not canonical.** The only invariant
    guaranteed is that `get_voigt`/`get_tensor` reconstruct the tensor. `dev` is *not*
    required to be trace-free, because two different pressure conventions coexist
    (`_trial_elastic_stress` uses `tr(ϵ)/3` even in 2D plane strain; the Voigt
    constructors use `-tr(σ)/S`, matching `σTr`/`yield_J2`'s divide-by-2-in-2D).
    `get_J2`/`get_τII` re-derive a trace-free deviator from `get_voigt` rather than
    trusting the stored `dev`. **Always consume a stored tensor via
    `get_voigt`/`get_tensor`, never by reading `.dev`/`.p` directly.**
  - **Return-mapping kernels stay in Voigt space internally** (`DP.jl`/`J2.jl` read
    `get_voigt(...)` once, run the closed-form/CPA algebra on `SVector`s, wrap once at
    the store) — deliberate, keeps numerics bit-identical to the pre-port code and
    preserves `elast_fast`'s hand-inlined performance.
  - Numerically verified bit-identical predictor math; the storage round-trip itself
    costs ≤1 ulp of the largest component. `test_performance.jl` measured -12.1%
    memory/-5.4% allocs/0.0% time vs. the pre-port baseline (this repo's typed-tensor
    objects allocate less than the bare `SMatrix`/`SVector` fields they replaced). JLD2
    needs no special handling — plain immutables over `T`/`SMatrix`.
  - `get_J2` (the tensor-invariant helper) was deliberately **not** unified with
    `retmap/J2.jl`'s same-named helper (which computes a different thing, `(‖ξ‖,n̂)`) —
    that one was renamed `_yield_normal` instead, so the collision is gone.
- ~~`PointSolidPhase.elast::E`/`AbstractElasticity`/`FiniteElasticity`/
  `LinearElasticity`~~ — **removed**, along with the `E` type parameter. `E` had
  decayed into a pure dispatch tag duplicating what `ST` already encodes once the
  tensor port landed (confirmed zero reads of `mpts.s.elast` before removal). Every
  kernel that dispatched on `E<:FiniteElasticity`/`E<:LinearElasticity` (`elast.jl`,
  `dynamic_relaxation/{fint,update}.jl`) now dispatches on `ST` instead.
- ~~`PointSolidPhase.P::Vector{T2}`~~ — **removed**, zero reads; pressure lives on the
  stress objects themselves (`mpts.s.σᵢⱼ[p].p`, or `get_voigt`'d).
- **`PointSolidPhase.ϵpII` is `Vector{SVector{2,T2}}`**, not `Matrix{T2}` — matches the
  rest of the codebase's per-particle `Vector{SVector}` convention. Since `SVector` is
  immutable, every in-place mutation site (`DP.jl`/`J2.jl`'s CPA updates, `nonlocal.jl`,
  `update.jl`'s bulk resets) reads the whole `SVector`, rebuilds it, writes it back —
  `drucker_prager`'s internal scratch stays a mutable `MVector{2,T}` for the actual
  return-mapping algebra; only the *stored* `mpts.s.ϵpII[p]` is immutable. Measured
  perf-neutral-to-better vs. the old `Matrix` layout (the `+=`-mutation pattern is
  ~25% *faster* as `Vector{SVector}`).
- **`get_voigt` (renamed from `get_vector`) is the single entry point for
  tensor↔Voigt-vector conversion, using the *engineering*-shear convention
  (`γxy=2εxy`)** — required wherever strain interacts with the elastic stiffness
  matrix `Del` (`σ=Del·ε_voigt` only works with the ordinary isotropic stiffness
  matrix under that convention). The old free functions `mutate`/`_mutate` existed
  purely to patch a factor-of-2 mismatch between `get_vector`'s old tensor-shear
  convention and `Del`'s engineering-shear expectation at each call site; both are
  gone now that `get_voigt` uses the correct convention directly.
  `LogarithmicStrain(ϵ::SVector)`/`InfinitesimalStrain(ϵ::SVector)` full-Voigt-vector
  constructors exist too, mirroring `CauchyStress(σ::SVector)`/`KirchhoffStress`;
  `get_voigt(LogarithmicStrain(ϵ)) == ϵ` to machine epsilon. `voigt_of(M::SMatrix)`
  (`elast.jl`) handles the one remaining case with no typed-tensor home (the Jaumann-
  rate correction term `σJ*ω'+σJ'*ω`, not itself "the" particle's stress).
  `retmap/MCRetMap.jl`'s own `mutate(...)` call is dead/unwired code, left calling an
  undefined function if ever exercised — not in scope to fix.
- **Unicode subscripts encode tensor rank — see "Conventions / house rules" below.**
  `PointSolidPhase.σᵢ`/`τᵢ` were renamed `σᵢⱼ`/`τᵢⱼ` for exactly this reason (a
  single-index subscript misleadingly read like a vector component on a full 2nd-order
  tensor field). **Gotcha discovered doing this rename**: a blind substring rename also
  caught genuinely-vector *local* variables (`get_voigt(...)` results inside
  `DP.jl`/`J2.jl`/`dynamic_relaxation/fint.jl`), which should have stayed single-index
  — fixed by renaming those specific locals back. Untangling this also surfaced a real,
  separately-introduced bug: `drucker_prager`'s return-mapped stress was briefly being
  written to a variable the function no longer returned, so it silently returned the
  pre-yield stress on every yielding step — verified fixed via a direct unit test (a
  hand-built past-yield stress state now returns an actually-changed result). Lesson:
  a mechanical rename like this needs a semantic check (tensor vs. vector), not just a
  substring match.
- **3D conformity check**: verify 3D simulations behave consistently across the
  explicit solver's configuration space (basis kind, `stab.locking`, `strain.deform`)
  — most verification to date has been 2D-only. 3D APIC currently hits a `DomainError`
  (sqrt of a negative value) even with `stab.locking=false`, root cause not chased
  (see "Thermal solution" below, which found and fixed unrelated 3D APIC bugs along
  the way but hit this separate, pre-existing instability).
- **`basis.how`/`basis.ghost` are GIMP-specific concepts living in the generic `basis`
  config section.** `how` is read unconditionally in `update.jl`'s dispatch regardless
  of basis kind; worth moving both onto `GimpBasis` itself (construction-time fields
  co-located with the kind that actually uses them) rather than basis-kind-agnostic
  top-level knobs `bsmpm`/`smpm`/`mlsmpm` silently ignore.
- **`Basis.N`/`∂N` storage rework**: currently plain `Matrix{T2}`/`Array{T2,3}`
  specifically to dodge the StaticArrays allocation-elision limits (see "Precision/
  StaticArrays gotchas" below) — correct and fast, but a workaround rather than a
  considered data-layout design.
- **`smpm`'s isolated gradient-consistency glitch**: `test_basis.jl` found
  `max|Σᵢ∂Nᵢ|` off by exactly `1.0` at one sweep point for `LinearBasis`, unlike
  `gimpm`'s broadly-degraded boundary PoU. Never root-caused — an error of exactly
  `1.0` smells like a real (if narrow) bug rather than float noise.
- **`DynamicRelaxationSolver`'s fate**: now that `solution` picks between
  `ExplicitSolver`/`ImplicitSolver` only, this third struct is orphaned scaffolding.
  Decide: wire it in as a third `solution` value, or state plainly that it's
  intentionally dormant.
- **Staged path to a poro-hydro-mechanical solution: mechanical solver (exists) →
  standalone fluid-only (hydro) solver → fuse the two.** Validate the hydro solver
  independently before coupling. `PointFluidPhase` (`lagrangian.jl`) needs real fields
  first — currently a literal empty placeholder, `mpts.f` always `nothing`.
- **Thermal solution — fixed and viable.** `thermal_problem` (renamed from
  `ic_thermal`) → `elastoplasm(jld2; workflows=[thermodynamic!])` runs end-to-end
  (verified: temperature stays bounded, cools monotonically from the boundary inward,
  both with `musl=true`/`false`). What was broken: stale `mpts.ϕ∂ϕ` reads across
  several kernels (replaced with the same `basis.N`/`∂Nrow(basis.∂N,...)` caching
  pattern the solid phase uses); a mistyped thermal `p2n!` overload dispatching on
  `mesh::Mesh` instead of `MeshThermalPhase`; no thermal-phase `mapsto` overload
  existed at all (added, mirroring the solid-phase one minus gravity/APIC/finite-
  strain concerns); `mesh.t`/`mpts.t` were unconditionally `nothing`
  (`setup_mesh`/`setup_mpts`/`setup_problem` gained a `thermal::Bool=false` kwarg);
  initial temperature was always zero regardless of config (off-by-one in
  `get_thermal.jl`). Not fixed, left as a known gap: `thermal_problem` defaults to
  `plot.status=false` — the plotting path has no `"T"` field entry and
  `display.jl` carries two incompatible `what_plot_field` dispatch styles, a separate
  piece of follow-up work. Separately found+fixed in the same pass: 3D APIC transfer
  for the *solid* phase had the identical dead-`ϕ∂ϕ` bug plus a mismatched mesh-phase
  type — fixed, but 3D+APIC then hits the pre-existing 3D instability noted under "3D
  conformity check" above; 1D APIC has a similar unfixed inconsistency (1D is
  essentially never exercised in this codebase).
- **Add a JLD2 time-series export option, as an alternative to plot-only checkpoints.**
  Every workflow's time loop calls `bake(mpts,mesh,t,solver)` at each checkpoint, but
  `bake` only renders a PNG when `solver.plot.status` is true — there's no way to
  record a field's evolution as raw data today. A natural fix reuses the same
  `checks` cadence and `solver.plot.what`-style field selection to append values into
  a JLD2 dataset at each checkpoint (`export.status`/`export.what`, sitting next to
  `bake()`). Writing the export is only half the feature — a companion loader/plotter
  for the saved time series would need to be built alongside it, not as a separate
  follow-up.

## Precision (32-bit) and StaticArrays performance gotchas

`dtype=(;T0=(Int32,Float32),bits=Int32(32),...)` is a supported but rarely-exercised
path — the default is always `(Int64,Float64)`, so bugs that only manifest under
`T1=Int32`/`T2=Float32` silently pass every day-to-day test. Two concrete patterns
found and fixed by actually running an `Int32`/`Float32` `slump_problem`→`elastoplasm`
end-to-end (not just skimming the code):

- **`@index(Global)` inside a `@kernel` function is always `Int64`**, regardless of the
  solver's index type `T1`. Passing it straight into anything dispatching on `::T1`
  throws a `MethodError` under `T1=Int32` — silently fine under the default `T1=Int64`
  only by coincidence. Always wrap it: `T1(p)`.
- **Bare numeric literals (`2.0`, `0.5`, `1/3`, ...) silently promote to `Float64`**,
  even inside a function generic over `T2`. Under `T2=Float32` it produces a mix of
  `Float32`/`Float64` values that eventually hits a `MethodError` several functions
  downstream of where the literal actually is. Always write `T2(2.0)` etc., never a
  bare literal, in any function generic over a float type parameter.

Separately, when caching per-particle results into a `Vector{SVector}`/`Vector{SMatrix}`
(or building one from scratch) inside a hot loop: **bundling `NN` computed values into
one whole-object `SVector{NN,...}`/`SMatrix{NN,D,...}` can silently heap-allocate**,
even when every intermediate step is individually zero-alloc in isolation — Julia's
escape analysis has a complexity/size cutoff past which it gives up eliding the
allocation, empirically true regardless of *how* the object was built (`MVector`
buffer, `ntuple`/`Val`, `hcat`, a flat tuple). Two mitigations that measurably worked
here (see `Basis.N`/`∂N`): (1) `ntuple(f, n)` needs `n` as `Val(n)`, not a plain `Int`,
to compile-time unroll — even when `n` is itself a type parameter, if used as a
runtime value it won't unroll; (2) store the *whole cached collection* as a plain
`Matrix`/`Array` written via scalar element assignment, rather than a `Vector` of
per-particle `SVector`/`SMatrix` objects — reading a row back out via `M[nn,:]` still
allocates for a runtime `nn`, so read it via a small `Val`-unrolled `ntuple` (see
`∂Nrow`) instead of array-slicing syntax. Always verify with `@allocated`/`@code_typed`
in a wrapping *function* (not top-level REPL scope, which has its own unrelated
allocation noise).

## Persistence (JLD2)

Simulation setup is saved/loaded as `ic["problem"]` (a `MechanicalProblem`, bundling
`mesh`/`mpts`/`time`) and `ic["basis"]` (its own top-level key, deliberately
independent of `Problem`); solver config as `cfg["solver"]` (NOT `cfg["instr"]` — that
key name is stale/wrong despite appearing in some scripts). `elastoplasm`/
`elastoplasm!` unpack `problem = file["ic/problem"]; mesh,mpts,time =
problem.mesh, problem.mpts, problem.time; basis = file["ic/basis"]`; `elastoplasm!`'s
write-back after each workflow rebuilds `file["ic/problem"]`/`file["ic/basis"]` in
place. Both take a `workflows::Vector{Function}` kwarg (plural — not `workflow`).
`mpts.s.cmp` (constitutive constants) persists as part of `ic["problem"].mpts`; there
is no separate `ic["cmpr"]` key. `test_column.jl` (formerly `test_collapse.jl`)
predated the current three-key layout for a long stretch (undefined `kwargser`, stale
`cfg["instr"]`, matrix-style `mpts.x` indexing) but was fully fixed — see Known bugs
history if picking that file up again ever surfaces a similar staleness.

## Operating the package

Typical end-to-end run (this is also the standard smoke test after any refactor):

```julia
using ElastoPlasm
L,nel = [64.1584, 64.1584/4.0], [40,10]      # domain size, element counts (2D here; 3 entries = 3D)
jld2  = slump_problem(L, nel; cli()...)           # builds mesh/mpts/basis/time, saves setup, returns the .jld2 path
out   = elastoplasm!(jld2; workflows = [elastodynamic!, elastoplastic!])
@assert out.success
```

Use `elastoplasm!` (not `elastoplasm`) whenever you need to inspect post-run
`mpts`/`mesh` state afterward — see Known bugs for why.

- `slump_problem(L, nel; fid="...", kwargs...)` (`src/home/script/example/slump.jl`) is
  the reference example problem. Calls, in order: `get_solver` → `setup_geometry` →
  `setup_problem` (internally: `setup_mesh` → `setup_material_constants` →
  `setup_mpts`, which calls `setup_cmp` to build `mpts.s.cmp` → `setup_time`; returns a
  `MechanicalProblem`) → `setup_basis`, then `export_problem(...)` to persist to
  `.jld2` and returns that path. Use it as the template for a new example problem.
- `collapse_problem(nel; kwargs...)` (`src/home/script/example/collapse.jl`) is a dry
  granular Drucker-Prager block-collapse example — a plain rectangular block released
  under gravity with no pre-existing slope, geometry via `get_collapse`
  (`src/home/init/mpts/get_collapse.jl`, masks the shared `mpts_populate` candidate
  grid down to a `[0,w0]×[0,h0]` box). Modeled directly on MaterialPointSolver.jl's
  `2d_collapse.jl` reference scenario (`LandslideSIM/Archive_MaterialPointSolver.jl_paper`)
  — domain/block size/material constants are hardcoded to match it, so `nel` is the
  only argument meant to normally change; everything else is an override. Uses the
  manual low-level `setup_geometry → setup_mesh → setup_material_constants →
  setup_mpts → setup_time → MechanicalProblem → setup_basis` pipeline (not the
  `setup_problem` wrapper), since it needs custom `ρ0`/`E`/`ν`/`ϕ`/`c0` that
  `setup_problem`'s fixed `(mesh,mat,solver)` geometry-helper signature can't pass
  through. Defaults `plast.status=true`/`plast.constitutive="DP"`, and — unlike the
  package default (`:roller`, frictionless normal-only slip) — a `:fixed` (sticky/
  no-slip) base boundary, needed for the material to actually pile up instead of
  sliding indefinitely; matches MaterialPointSolver.jl's own boundary treatment.
  **Not to be confused with `column_problem`** (`src/home/script/example/column.jl`,
  renamed from the old `collapse_problem`) — a 1-D elastic self-weight column
  convergence test with no plasticity at all; see the Known bugs entry on
  `test_collapse.jl`/`test_column.jl` for that history.
- `cli()` (`src/home/api/solver/cli.jl`) parses interactive/CLI overrides for solver
  config; `cli(; ui=true)` prompts interactively, `cli()` with no args picks defaults
  non-interactively — pass its result as `kwargs...` into `slump_problem`/`get_solver`.
- `get_solver(; dim=2, kwargs...)` (`src/home/api/solver/get_solver.jl`) merges
  user kwargs over `get_default()`, resolves the execution backend, and builds the
  `Cairn` dispatch table (`cairn.ignite`/`.mapsto`/`.update`/`.implicit`).
- `elastoplasm(jld2_path; workflows=[...])` reopens the `.jld2`, unpacks
  `(mesh, mpts, basis, time, solver)`, and runs each workflow function in order as
  `workflow!(mpts, mesh, basis, time, solver)`. `elastoplasm!` is the in-place variant
  (opens `"r+"` so postprocessing can write back). Both accept any callable matching
  that signature — the built-ins are `elastodynamic!`/`elastoplastic!` (explicit
  solver) and `elastoquasistatic!` (dynamic_relaxation), freely composable in the
  `workflows` vector since they share the same call signature.
- Plotting is opt-in via `solver.plot.status`; when on, `elastoplasm` auto-saves a PNG
  per workflow named `"$(dim)_$(basis.which)_$(solvertype)_$(deform)_$(workflow)_$(quantity).png"`
  under the run's `dump/.../plot` path.

### Running tests / benchmarks

```julia
julia --project=. -e 'using Pkg; Pkg.test()'          # or: include("test/runtests.jl")
```

`test/runtests.jl` interactively lets you pick which `test/testset/test_*.jl` files to
run (auto-runs all of them under `GITHUB_ACTIONS=true`). It swallows exceptions inside
`runtests()`'s `try/catch` (just increments a fail counter) — when debugging a failing
test file, `include()` it directly instead of through `runtests()` to see the real
stacktrace, e.g.:

```julia
using Test,JLD2,ProgressMeter,Suppressor,Plots,LaTeXStrings,REPL.TerminalMenus,ElastoPlasm
using BenchmarkTools, KernelAbstractions
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync
global ROOT, DATASET = joinpath(pwd(),"test"), joinpath(pwd(),"test","dataset")
include("test/testset/test_performance.jl")
```

`test_performance.jl` benchmarks core kernels directly via `solver.cairn.*` (including
`ignite.shpfun!`) and compares against `test/dataset/performance_baseline.jld2` (keyed
by CPU name); delete/update that file deliberately if a change is a genuine, intended
perf shift rather than a regression.

`test_basis.jl` builds a small `n×m` mesh and, for every basis kind, sweeps a particle
along the mid-height node row, plotting each tracked node's `Nᵢ(x)`/`∂Nᵢ/∂x(x)` and the
partition-of-unity `Σᵢ Nᵢ(x)` — useful as a quick correctness check on a new/changed
basis kind. Both this file and `test_performance.jl` report numeric checks
informationally (`println`, not hard `@test` gates), deliberately — not every basis
kind is expected to hold exact PoU everywhere, so a strict `@test` would be a false
failure rather than a real regression signal. When adding LaTeX to a generated figure:
wrap genuine math expressions but keep plain identifiers (basis kind names) as plain
strings — a multi-letter plain word in math mode doesn't render reliably under the
default GR backend.

### Controlling solver behaviour via `defaults.jl`

`get_default()` (`src/home/api/solver/defaults.jl`) returns the base `NamedTuple`
config merged by `get_solver` — this is the single place that defines every tunable
knob and its out-of-the-box value.

- `solution` — `"explicit"`/`"implicit"`, picks the concrete solver struct at
  construction time (`ExplicitSolver`/`ImplicitSolver` — same field layout,
  `DynamicRelaxationSolver` unused scaffolding). Also what `elastoplasm.jl`'s
  generated filenames and `logs.jl`'s startup banner print. Picking `solution=
  "explicit"` then calling an implicit-only workflow fails with `MethodError` at the
  workflow call, not silently — `elastodynamic!`/`elastoplastic!` are typed to
  `ExplicitSolver` specifically. `ignite()` is the one place dispatching on the
  abstract `AbstractSolver` instead, since both paths call it regardless of `solution`.
- `dtype` — arithmetic precision (`bits`, element types `T0`)
- `basis` — `which` (`"bsmpm"`/`"gimpm"`/`"smpm"`/`"mlsmpm"`, see `get_basis`), `how`
  (GIMP domain update mode), `trsfr` (P2G/G2P transfer scheme: `"std"`/`"tpic"`/
  `"apic"`, see `get_transfer`), `C_pf` (PIC/FLIP blend). `trsfr`/`C_pf` live here
  rather than a separate `transfer` section since `Basis` (the struct) owns both
  `kind` and `transfer` as sibling fields — see "Transfer scheme dispatch" above.
- `strain` — `deform` (`"finite"`/`"infinitesimal"`)
- `stab` — `locking` (F-bar volumetric locking correction on/off), `damping`, `musl`
  (MUSL velocity reprojection on/off — lives here rather than under `basis`/`transfer`
  since it's a stabilization technique applied regardless of transfer scheme).
- `bcs` — `dirichlet` boundary condition matrix, one `[lower upper]` row per dimension
- `grf` — Gaussian random field generator for heterogeneous cohesion/friction fields
  (`status` toggles it on; see `GRF.jl`)
- `plast` — `status`, `constitutive` (`"DP"`/`"VM"`/`"MC"`/`"camC"` — not all are
  wired up, check `setup_cmp`'s branch before relying on one; the `retmap` kernel
  dispatch lives on `Point`'s `CM` type parameter, not a runtime string — see "DP/J2
  retmap kernel unification")
- `nonloc` — non-local plastic strain regularization (`status`, `ls` length scale)
- `plot` — `status`, `freq` (plot every N `Time` checkpoints), `dpi`, `what` (list of
  field specs to plot, keyed by name via `get_mpts_variable_config()`)
- `perf` — `status`; when true, forces `deform="infinitesimal"` and disables `nonloc`
  and swaps in the `_fast` kernel variants for a lighter-weight run
- `backend` — `select` (execution backend, `"host"`/GPU target) and `distributed`

Two ways to change solver behaviour:

1. **Per-run override (preferred, no code change)** — pass any of these keys as
   keyword arguments through `slump_problem`/`get_solver`; `get_solver` only keeps
   kwargs whose keys already exist in `default` (unrecognized kwargs get a `@warn`,
   not an error), then merges them over `default` via `merge(ref, user)`. **This is a
   shallow-per-key merge** — to override one field of a nested `NamedTuple` block
   (e.g. just `basis.trsfr`) you must pass the *whole* section
   (`basis = (; which=..., how=..., trsfr=..., C_pf=...)`), not just the one field,
   since `merge` replaces the whole `:basis` entry rather than recursing into it —
   this has caused real, previously-live bugs (`thermal_problem`'s and
   `test_basis.jl`'s own hardcoded `basis=(;which=...)` overrides silently dropped
   `trsfr`/`C_pf` until fixed to spread `get_default().basis...` first). Example:
   ```julia
   slump_problem(L, nel; basis=(;which="gimpm",how="Uii",trsfr="apic",C_pf=1.0), strain=(;deform="finite"), stab=(;locking=true,damping=0.1,musl=true))
   ```
2. **Change the package-wide default** — edit the literal value in `get_default()`
   directly. Do this only for a genuine change of the shipped default behaviour, not
   for one-off experimentation (use kwargs for that).

New solver-config knobs must be added as a new key inside the relevant block (or a new
top-level block) in `get_default()`'s `default` NamedTuple, then threaded through
wherever `init_ignite`/`init_mapsto`/`init_update`/`init_implicit` branch on
`instr[:section][:key]` to select kernels — grep those `init_*` functions for the
existing pattern before adding a new one. If the new field belongs on the solver
struct itself (like `solution` did), it also needs threading through
`ExplicitSolver`'s/`ImplicitSolver`'s field list (`src/boot/needs/types/solver.jl`) and
the corresponding positional arg in `get_solver.jl`'s final
constructor call — both structs are kept in lockstep by hand.

### Using `cli()`

`cli(; ui::Bool=false)` turns `get_default()`'s config tree into a `Dict{Any,Any}` of
kwargs suitable for splatting into `slump_problem`/`get_solver` as `cli()...`.

- `cli()` (default, `ui=false`) — **non-interactive**, returns `get_default()`'s values
  reshaped into the kwargs `Dict`. What you want for scripts/automated runs.
- `cli(ui=true)` — **interactive**: first a `MultiSelectMenu` to pick which top-level
  sections to configure, then a `RadioMenu`/`MultiSelectMenu` per leaf option (nested
  sections handled recursively; sections with a `status` field skip sub-prompts when
  answered `false`; `plot.what` gets multi-selection). Unselected sections fall back
  to `get_default()`. Use only from a real terminal/REPL — blocks on
  `TerminalMenus.request`, never call from a non-interactive script or agent session.
- To override specific options while using `cli()` for everything else, splat `cli()`
  first then override afterward — later kwargs win. Same shallow-merge caveat applies.
- `get_option()` is a reference for valid values per key (e.g.
  `get_option().basis.trsfr` → `["std", "tpic", "apic"]`) — some tunables (e.g.
  `plast.constitutive`) accept values it doesn't fully enumerate; check `init_update`'s
  dispatch in `update/update.jl` if in doubt whether a value is actually wired up.

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
  or — as happened this session — silently operating on the untouched pre-simulation
  state because `elastoplasm` was used instead of `elastoplasm!`).
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
  entry under Planned improvements for a concrete example: local Voigt-vector
  variables getting swept up in a tensor-field rename).

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
