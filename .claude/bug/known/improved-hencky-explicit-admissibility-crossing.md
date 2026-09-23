# `material.elastic="improved_hencky"` crashes the explicit slump with plasticity

**Status: open.**

`slump_problem` (2D, `[40,10]`, `cli()` defaults, `material.plastic="DP"`) run through
`elastoplasm!(...; workflows=[elastodynamic!, elastoplastic!])` fails under
`"improved_hencky"`. The `"hencky"` and `"hypoelastic"` runs of the same setup finish cleanly.
`elastodynamic!` alone runs fine under `"improved_hencky"`; the failure only occurs in the
`elastoplastic!` stage.

What is known:

- **Timestep fixed, but that alone doesn't solve it.** `get_dt` used the constant `Kc` for the
  wave speed, while improved Hencky's tangent bulk modulus is `Kc/n` at rest (10× `Kc` at the
  default `n₀=0.1`) and diverges as `n→0⁺` (Pretti et al. 2024, Eqs. 36-37). `get_dt` now
  dispatches `Ktan` on `Point`'s `SM` and uses Eq. 36 for `ImprovedHenckySolid`. Before this,
  the run died in F-bar (`ΔJ<0`, `volumetric.jl`). Lowering the CFL factor to 0.05 by hand
  also avoided it. After the fix it still dies, now in `get_dt` itself: `Ktan<0`.
- **Root cause: particles cross the admissibility bound.** `Ktan<0` means `n<0`, i.e. a
  particle's elastic `ϵᵥ` fell below `ln(1-n₀)` (paper's Eq. 30a, about 10.5% volumetric
  compression at `n₀=0.1`). An explicit step jumps across the asymptote the paper relies on
  (it uses an implicit solver), after which `n`, the pressure and `Ktan` all flip sign.
- **The slump itself compacts that much.** The `"hencky"` run of the same slump ends with
  total porosity `mpts.n` clipped to `0` (`deform.jl`) for some particles, so this is not
  specific to the improved law.

Open decision, not attempted: whether the fix is
- a physically realistic `n₀` for the slump,
- an explicit error when Eq. 30a is violated (no silent clamp),
- or a different step control.

Do not clamp `n` silently — that hides the violation.
