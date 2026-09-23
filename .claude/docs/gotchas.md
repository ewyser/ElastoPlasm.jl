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
here (see `Basis.N`/`∂N` in `architecture.md`): (1) `ntuple(f, n)` needs `n` as
`Val(n)`, not a plain `Int`, to compile-time unroll — even when `n` is itself a type
parameter, if used as a runtime value it won't unroll; (2) store the *whole cached
collection* as a plain `Matrix`/`Array` written via scalar element assignment, rather
than a `Vector` of per-particle `SVector`/`SMatrix` objects — reading a row back out
via `M[nn,:]` still allocates for a runtime `nn`, so read it via a small `Val`-unrolled
`ntuple` (see `∂Nrow`) instead of array-slicing syntax. Always verify with
`@allocated`/`@code_typed` in a wrapping *function* (not top-level REPL scope, which
has its own unrelated allocation noise).
