# Drucker-Prager tension/apex return produced `(σm−P)·I` and flipped the deviator beyond the apex

**Status: fixed.**

`_druckerprager_return_map` (`retmap/DP.jl`) handled trial states at or past the tension
cutoff in two ways, both wrong for this code's configuration:

- **Case (d), `h≤0 && P≥σm`:** it set `Pn = σm - P` and passed it to `σn`, which adds `Pn` to
  the diagonal as an *absolute* mean stress. The result was `σ = (σm−P)·I`, i.e. a compressive
  mean for any tensile overshoot, and discontinuous with the cone branch at `P = σm` (cone gives
  `σm`, apex gave `0`). It looks like Huang et al. 2015 Eq. 55's *increment* `(σt − σm*)` was
  ported as the absolute value.
- **Case (e), `h>0 && P≥σm`:** it went through the smooth-cone formulas. With the hardcoded
  dilatancy `ψ=0` (`ηB=0`) those keep the mean at `P ≥ σm` and give `τn = ξ − η·P ≤ 0`, so `σn`
  scaled the deviator by a negative factor, flipping its sign.

Root cause: the code places the tension cutoff exactly at the apex (`σt = σm = ξ/η`, Huang's
Eq. 19 with `σt = σt_max`, which Huang himself uses when `kφ = 0`). In that degenerate case,
neither of Huang's single-surface returns is admissible: literal Eq. 55 keeps `τ*` at the apex
(so `fs = τ* > 0`), and case (e) flips the deviator as above. `agent-mpm-specialist` re-derived
both from Huang's Eqs. 17-24, 43, 45-47 and 53-56 and confirmed them.

Fix: for `P ≥ σm`, the apex `σm·I` is the only admissible state (`fs ≤ 0` forces `τ = 0`,
`ft ≤ 0` forces mean `≤ σm`), so both cases now return there. `Δλ = (P−σm)/Kc` (Eq. 54) and
the `ϵpII` increment `√2·Δλ/3` (Eq. 56) are unchanged, and the cone branch is gated on
`P < σm` only. `σn(σm, τ0, 0, 1)` is used so a purely hydrostatic trial state (`τII = 0`)
no longer gives `0/0 = NaN`.

Verified:
- **Direct unit calls on hand-built trial states** (cases d and e, hydrostatic, `σm±1 Pa`):
  every return lands exactly on `σm·I` with `fs = ft = 0`, stays continuous across `P = σm`,
  gives no NaN, and leaves elastic states untouched.
- **Slump with DP (2D `[40,10]`):** `hypoelastic` is bit-identical, since it never reaches the
  apex. `hencky` changes slightly on every particle (`max|u|` 0.01226 → 0.01197, `ϵpII` max and
  mean unchanged).

Deferred, not fixed: at the apex the `ϵpII` increment carries only the tensile part. A
two-surface corner return would add the shear-flow part `√(1/3)·Δλˢ`, `Δλˢ = τII/(2Gc)`, which
makes `Δλ` jump from ~1e-2 just below the apex to ~1e-7 just above it. That rests on
multi-surface (Koiter) corner-return theory, e.g. de Souza Neto, Perić & Owen 2008 §8.3, which
is not in `refs/`, so it waits until that source is added. `Borja_2013.pdf` was checked; it only
covers π-plane corners.
