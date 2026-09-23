# `material.elastic="improved_hencky"` crashes the explicit slump with plasticity

**Status: open.**

`slump_problem` (2D, `[40,10]`, `cli()` defaults, `material.plastic="DP"`) run through
`elastoplasm!(...; workflows=[elastodynamic!, elastoplastic!])` fails under
`"improved_hencky"`. The `"hencky"` and `"hypoelastic"` runs of the same setup finish cleanly.
`elastodynamic!` alone runs fine under `"improved_hencky"`; the failure only occurs in the
`elastoplastic!` stage.

What is known:

- **Timestep: fixed, necessary but not sufficient.** `get_dt` used the constant `Kc` for the
  wave speed, while improved Hencky's tangent bulk modulus is `Kc/n` at rest (10× `Kc` at the
  default `n₀=0.1`) and diverges as `n→0⁺` (Pretti et al. 2024, Eqs. 36-37). `get_dt`
  (`explicit/get.jl`) now calls `_Ktan(EL, strain, cmp, n₀)`, dispatched on `Point`'s `EL`,
  with Eq. 36 for `ImprovedHenckySolid` and `Kc` otherwise. Before this, the run died in F-bar
  (`ΔJ<0`, `volumetric.jl:29`); lowering the CFL factor to 0.05 by hand also avoided it. After
  the fix it still dies, now in `get_dt` itself: `_Ktan<0`, i.e. `n<0`, i.e. a particle's
  elastic `ϵᵥᵉ` fell below `ln(1-n₀)` (paper's Eq. 30a, ≈-0.105 at `n₀=0.1`).
- **Likely root cause (identified by reading, not yet confirmed by a fix): the return map
  recovers the elastic strain with the *linear*-Hencky inverse.** Both `retmap/DP.jl:97` and
  `retmap/J2.jl:83` rebuild the post-yield strain as `LogarithmicStrain(cmp.Del\τᵢ)`, i.e.
  `ϵᵥᵉ = -P/Kc`. Under improved Hencky the actual relation is `P ≈ Kc·ϵᵥᵉ/n` (Eq. 34), so the
  recovered `ϵᵥᵉ` is too large in magnitude by roughly `1/n` (≈10×). In compression, a single
  yielding step throws `ϵᵥᵉ` past `ln(1-n₀)`; on the next step
  `n = 1-(1-n₀)/exp(ϵᵥᵉ)` goes negative and the pressure and `_Ktan` flip sign. This matches
  every observation: elastic-only runs survive, plastic runs don't, and a smaller CFL only
  delays it. The `#= WIP ... =#` block right below the `DP.jl` line already sketches a
  tangent-based alternative.
- **Proposed fix, grounded in the paper rather than improvised:**
  - **Smooth-cone return:** `DP.jl` runs with dilatancy `ψ=0`, so `ηB=0` and `Pn=P`. The
    volumetric elastic strain must stay at its trial value; only the deviator changes,
    `dev = τ_dev/(2Gc)`, which Eq. 33 leaves linear.
  - **Apex return:** the pressure changes, so `ϵᵥᵉ` has to come from inverting Eq. 34, which is
    nonlinear. Transcribe the paper's Appendix A.2 (Eq. A.9 system) for this rather than
    inventing a scheme.
  - Implement the smooth-cone part first and check whether the slump runs before touching the
    apex branch.
- **Superseded explanation.** This file first blamed the slump's own compaction crossing
  Eq. 30a. The `"hencky"` run of the same slump does end with total porosity `mpts.n` clipped
  to `0` in `deform.jl` for some particles, so the slump compacts strongly either way. That
  alone doesn't explain a crash that only happens with plasticity; the strain recovery above
  does. Keep in mind while fixing: do not clamp `n` silently, since that hides a real Eq. 30a
  violation.
