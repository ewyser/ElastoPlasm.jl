# Non-local plastic-strain regularization (`nonlocal.jl`) was `O(nmp²)` in both memory and time, every timestep

**Status: fixed, and a real correctness bug found along the way.**

`Basis.e2p`/`p2p` (`nmp×nel`/`nmp×nmp`) were **unconditionally dense**
(`setup_basis.jl` built them via `T1.(spzeros(...))`, but the struct fields
were plain `Matrix`, so the sparse-zeros broadcast got densified on
construction) regardless of `solver.nonloc.status` — same class of bug as
`bug/fixed/mesh-solid-phase-dense-mij-matrix.md`, and just as live:
`nonloc.status` **defaults to `true`** (`defaults.jl`), independent of
`plast.status`, so this ran by default in essentially every simulation.
Every timestep also allocated a fresh `w = spzeros(T2,nmp,nmp)` and wrote
into it via scalar `w[p,q]=...` — `setindex!` on a `SparseMatrixCSC` at a
previously-zero location is `O(nnz)` (shifts every higher-index stored
nonzero), making this pattern worse than even a dense fill. Measured at
`nmp=10133`: 2.213 GiB → 986 MiB per `elastodynamic!` run (`plast.status=
true`, 1137 steps); at `nmp=40413` the fixed version completed in
97.6s/1.433 GiB while the pre-fix version was still running after 30+
minutes (not fully timed to completion — analytically it's the expected
`O(nmp²)` inner-loop-iteration blowup: the old `"p->q"`/`"p<-q"` phases
scanned **all** `nmp` particles per yielding particle instead of using
`basis.e2e` — already a genuinely sparse, `ls`-bounded element-neighbor
structure — to narrow the search).

Fix: a per-step CSR element→particle bucket list (`build_el2p`, `tplgy.jl`,
`O(nmp+nel)`, sequential) combined with `basis.e2e` bounds the neighbor
search to `O(nmp×k)` (`k` = local particle density within `ls`,
resolution-independent); `Basis.e2p`/`p2p` and the per-step `w` matrix are
gone entirely (pairwise weights are cheap enough — one `sqrt`+`exp` — to
recompute in both passes rather than cache). The old single kernel branching
on a `type::String` ("tplgy"/"p->q"/"p<-q") was also split into two
separately-named kernels, `nonlocal_pq`/`nonlocal_qp`
(`update[:nonloc_pq!]`/`[:nonloc_qp!]` in the `Cairn` table) — matching the
one-kernel-per-concern convention used everywhere else in this codebase
(`elast!`, `ignite`/`mapsto`'s per-phase entries) and dropping a per-particle
`String` comparison that no longer served any purpose once the `"tplgy"`
phase was gone.

## Correctness bug found while verifying the fix

**While verifying this fix bit-for-bit against the old code, found the old
code was never actually correct**: it deduplicates each unordered pair
`(p,q)` via `if w[p,q]==0` before recording it into `p2p`, processing
particles in ascending index order — whichever particle in a pair is
visited *first* claims the pair and writes `p2p[q,p]`; the *other*
particle's later iteration finds `w` already nonzero and never writes its
own side. Net effect: for a mutual-neighbor pair, only the *lower-indexed*
particle ends up receiving that neighbor's contribution in the `"p<-q"`
averaging pass — the *higher-indexed* particle silently loses it. Confirmed
directly: with every particle yielding, 8086 of 8677 recorded neighbor-pair
entries in the old `p2p` were asymmetric (present in only one direction) —
the regularization has been missing roughly half of each particle's true
neighbors, systematically, since this code was written. The rewrite computes
each particle's own neighbor sum independently (no cross-particle writes, no
dedup guard) and is verified symmetric and correct by independent
brute-force recomputation of `W`/`ϵpII[2]` for several particles — it does
**not**, and structurally cannot, reproduce the old asymmetric numbers; that
is the fix, not a regression.
