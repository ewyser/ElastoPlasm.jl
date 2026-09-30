# Hypoelastic `elast!`: the Jaumann rotation term is identically zero

**Status: fixed** (branch `fix-hypoelastic-jaumann-rotation`).

The small-strain (`material.elastic = "hypoelastic"`) kernel in
`src/home/core/solver/explicit/update/fun/elast.jl` adds the objective-rate correction as

```julia
jaumann = CauchyStress(σJ * ω' + σJ' * ω)
```

with `σJ` the (symmetric) Cauchy stress and `ω = ½(∇v − ∇vᵀ)` (antisymmetric). Since `σJ' = σJ`
and `ω' = −ω`, this is `−σω + σω = 0` for every input — checked numerically (a random symmetric
`σ` and antisymmetric `ω` give exactly `0.0`, while `ωσ − σω` is O(1)). Hypoelastic runs therefore
integrate `σ ← σ + Del·ϵ` with no rotation at all: rigid-body spin rotates nothing.

**History.** The generic kernel carried the same expression before `c9da3a7b` ("Massive changes in
implementation logic", 2026-08-09), in matrix form (`σJ*ω' + σJ'*ω`), and `dev` still has it. The
hand-unrolled `_fast` kernels did rotate the stress: in 2D `σxx += 2ω_xy σxy`,
`σyy −= 2ω_xy σxy`, `σxy += ω_xy (σyy − σxx)`, i.e. `ωσ − σω` component by component. Those were
only reachable through `perf.status`, which also forced `hypoelastic`; when they were removed (see
"`_fast` kernels and the `perf` config section — removed" in `planned-improvements.md`), the
equivalence check ran on the default Hencky path, so the loss of the rotation in hypoelastic runs
went unnoticed.

**Fix.** `jaumann = Δω·σ − σ·Δω` with the spin increment `Δω = ½(ΔF − ΔFᵀ)`. The rate form is
Borja (2013), Eq 5.46 (`τ∇ = τ̇ − w·τ + τ·w`, so `τ̇ = τ∇ + wτ − τw`, `w = skew(∂v/∂x)`); the
incremental form matches Huang et al. (2015), Eqs 29/32 (`σᴿ = σ + σᵢₖωⱼₖ + σⱼₖωᵢₖ`), both in
`refs/`. `∇vᵢⱼ = ∂vᵢ/∂xⱼ` is a rate and `ΔF = I + Δt·∇v` (`deform.jl`), so `½(ΔF − ΔFᵀ) = Δt·w`
is the increment `elast!` needs without taking `Δt`, from the same `ΔF` as the strain increment.
It stays exact on the in-plane D×D block: `Δω` has no out-of-plane entries and the 2D stress has
no out-of-plane shear. The stored `mpts.s.ωᵢⱼ` (the rate, written but never read) is unchanged.

**Verification.** Rigid rotation through the `elast!` kernel: a pre-stressed plane-strain state
(`σzz` included) rotated by 90° in N steps of `ΔF = R(Δθ)`. Before, `|σ − Rσ₀Rᵀ|/|σ₀|` saturated at
1.13 as `Δθ → 0` (the stress never rotated: `σxy` stayed +25 kPa, exact −25 kPa); after, it falls
to 1.5e-2 at `Δθ = 5e-4`, first order in `Δθ` (the explicit update plus the `cos Δθ − 1` strain an
exact rotation `ΔF` carries). Hypoelastic 2D slumps, before → after: elastic max |Δx| 2.0 cm,
von Mises 16 cm, Drucker–Prager 59 cm (medians 0.4, 2.2, 3.4 cm); stress ≤4.4% of its maximum; all
runs complete. Hencky runs are untouched. Not run: 3D hypoelastic (the same formula on 3×3) and the
`test_workflow` matrix.

**Consequence for the storage refactor.** On `feat-refactor-tensor-storage` the term is built on the
in-plane block through `CauchyStress(M::SMatrix{D,D})`; results were verified bit-identical, so
the refactor preserved the bug, it did not introduce it; the fix was made on top of it.
