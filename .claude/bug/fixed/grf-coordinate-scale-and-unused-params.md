# Gaussian random field used halved coordinates, an unconverted correlation length, and ignored `grf.param`

**Status: fixed.**

The cohesion random field (`GRFS_gauss`, `init/mpts/GRF.jl`) had three defects relative to its
source, Räss, Kolyukhin & Minakov (2019), *Computers & Geosciences* 131:158-169
(`refs/ElastoPlasm/grf/Raess_etal_2019_random3D.pdf`):

1. **Halved coordinates.** The phase used grid step `dx = 0.5·Δx/ni`, while particles sit `Δx/ni`
   apart, so the field was evaluated at half the physical coordinates. `dy`/`dz` also used the
   x-spacing `Δx`. The paper's Algorithm 3 evaluates the phase at physical coordinates. The
   maintainer confirmed the ½ was unintended.
2. **Correlation length not converted.** Wave numbers were divided by `I` directly. Eq. (13)
   defines `C(r) = σ²exp(-r²/l_f²)` with correlation length `I = l_f√π/2`, so Algorithm 2 uses
   `k·√2/l_f` with `l_f = 2I/√π`.
3. **`grf.param` ignored.** `Iₓ`, `Nₕ`, `kₘ` (offered by `get_default`/`cli`) were never read;
   `sf = 5e3`, `If = 2.5`, `Nh = 5000`, `k_m = 100` were hard-coded.

Fix: `grf_gauss(xp, σ, I, Nh, kₘ)` (since renamed `random_field`, which also covers the exponential covariance) implements Algorithm 2 at the material points' own coordinates
(D×N), called by `get_cohesion(xp, mat, solver)` on the particles each problem actually keeps
(`get_slump`, `get_collision`), instead of on the candidate grid. `mpts_populate` now returns
positions only. `grf.param` gains `σ` (default `5e3`, the old hard-coded value) and all
parameters are read from it. `"exponential"` is no longer offered by `get_option` and errors
clearly if set. Per-axis `Iₓ` scale each wave-vector component by its own `l_f`: that's an
extension of the paper's isotropic Gaussian case, identical to it when all three are equal.

Verified:
- **Against Eq. (13):** on an 80×80 m point cloud, the empirical correlation at r = 1, 2, 3 m is
  0.885/0.599/0.319 against the analytic 0.882/0.605/0.323 (`I = 2.5`), in both x and z. An
  anisotropic `I = [5,5,1]` follows each axis's own `l_f`.
- **Snapshot:** every example without the random field is bit-identical; the random-field slump
  changes only `mpts.s.cmp`, and runs end to end.
- **Speed:** a type-stable phase loop gives 26k points × 5000 harmonics in 2.4 s.
