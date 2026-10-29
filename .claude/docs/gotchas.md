## Precision (32-bit) and StaticArrays performance gotchas

`dtype=(;T0=(Int32,Float32),bits=Int32(32),...)` is a supported but rarely-exercised
path — the default is always `(Int64,Float64)`, so bugs that only manifest under
`T1=Int32`/`T2=Float32` silently pass every day-to-day test. Two concrete patterns
found and fixed by actually running an `Int32`/`Float32` `slump_case`→`execute!`
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
here (see `Basis.N`/`∂N` in `architecture.md`): (1) `ntuple(f, n)` needs `n` as
`Val(n)`, not a plain `Int`, to compile-time unroll — even when `n` is itself a type
parameter, if used as a runtime value it won't unroll; (2) store the *whole cached
collection* as a plain `Matrix`/`Array` written via scalar element assignment, rather
than a `Vector` of per-particle `SVector`/`SMatrix` objects — reading a row back out
via `M[nn,:]` still allocates for a runtime `nn`, so read it via a small `Val`-unrolled
`ntuple` (see `∂Nrow`) instead of array-slicing syntax. Always verify with
`@allocated`/`@code_typed` in a wrapping *function* (not top-level REPL scope, which
has its own unrelated allocation noise).

## Thread count and bit-identity checks

Multi-threaded runs are not reproducible to the last bit: the atomic nodal sums of P2G add in
thread order. With 8 threads the §4.3 block under dynamic relaxation (thousands of iterations)
differs by ~1e-10 m in position between two runs of the *same* code, while the explicit slump
happened to agree. So compare a refactor against its parent with `-t 1`: single-threaded, the
block matched the pre-refactor commit exactly (`refactor-init-problem-spec`). Setup-only checks
(`test_setup.jl`) are not affected.

**GIMP's stretch domain update is sensitive to how `deform` compiles.** `particle_domain` (`common/domain.jl`)
computes `eigen(Symmetric(F'F))` and `n*diagm(√λ)*n'` with StaticArrays, whose products use `muladd`: LLVM fuses
them into FMAs or not depending on the code around them, so a change elsewhere in the `deform` kernel can move `ℓ` by
one ulp per step, and the whole GIMP (stretch) run by round-off. Seen twice: `fix-gimpm-stability` (moving the update
into `deform`), and `feat-gpu-explicit`'s parametrized containers (`MeshProperties.nel` a static vector, so the
kernel's loop bounds became constants): the first difference was `mpts.ℓ` at step 2, nothing upstream of it, and the
GIMP·APIC·von Mises·hypoelastic slump ended at most 5e-14 relative apart (stress), 2e-16 (positions); uGIMP, which has
no update, stayed bit-identical. In a bit-identity check, treat this case's difference as round-off once its first
divergence is `ℓ`.


## Pkg on this machine rewrites the Manifest's Julia version

The Manifest is written by Julia 1.10.10; this machine runs 1.10.5. Any `Pkg.develop`/`Pkg.resolve` here sets
`julia_version = "1.10.5"` and downgrades the stdlib `OpenLibm_jll` (0.8.5+0 → 0.8.1+2). Restore those two lines
before committing, so that the Manifest diff is only the intended change (adding EasyDraw by path was).

## Working with the submodules

`lib/EasyDraw` and `lib/UnifiedBackend` are git submodules: checkouts of EasyDraw.jl's and UnifiedBackend.jl's own
repositories, of which ElastoPlasm records one commit each. Below, `lib/X` is either.
- **Clone** with `git clone --recursive`, or run `git submodule update --init` after a plain clone; without it
  `lib/X` is empty and the package does not load (the Manifest finds both by path).
- **Change one:** in `lib/X`, switch to its working branch first (a submodule is checked out on a detached commit;
  EasyDraw's is `main`, UnifiedBackend's `dev`, merged into its `main` by pull request), commit and push there, then in
  ElastoPlasm `git add lib/X` and commit, which records the new commit. Push the submodule before ElastoPlasm, or others
  get a commit they cannot fetch. `git status` in ElastoPlasm shows `modified: lib/X (new commits)` until then.
- **Take another commit:** `git -C lib/X pull` (or `git submodule update --remote lib/X`), then `git add lib/X` and
  commit; after switching ElastoPlasm branches, `git submodule update` checks out the commits that branch records.
- **Branches from before a submodule** (`dev`, …) do not have it: switching to one leaves the folder there, untracked;
  it does no harm and comes back into use when switching back.
- A sibling clone of UnifiedBackend.jl (`../UnifiedBackend.jl`) is a separate checkout: `git pull` it after pushing
  from `lib/UnifiedBackend`.

## Device runs differ from host runs by round-off

StaticArrays' small matrix products use `muladd`, which LLVM fuses into FMAs or not depending on how the kernel is
compiled. A kernel compiled for a device's array types fuses differently from the host's: on JLArrays (the local
stand-in device, `test_device.jl`), `F = ΔF·F` in `deform` first differed by one ulp while `∇v` and `ΔF` were equal,
and a short plastic slump ended at most 9e-12 relative apart (Float64), 8e-7 (Float32); integer and Bool arrays were
equal. Compare device and host runs with a relative tolerance, not bit for bit; a GPU's own `exp`/`log`/`sqrt` add to
it. JLArrays has no `KernelAbstractions.synchronize` method (the test defines a no-op one).
