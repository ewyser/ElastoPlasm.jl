# Von Mises return mapping was not radial (shear flow direction off by a factor of 2)

**Status: fixed** on `feat-implicit-dr-sutcliffe` (`fix: Make the von Mises return mapping radial`).

`_vonmises_return_map` (`retmap/J2.jl`) iterated `Δσ = Δλ·Del·n` and `Δλ = f/(nᵀ·Del·n)` with
`n = ξ/‖ξ‖`, the unit deviator in *stress*-Voigt form (shear stored once). `Del` acts on
*engineering*-Voigt strains, so the shear components of the flow direction must be doubled.
Without that, normal components were corrected with `2G` and shear with `G`: the CPA still
converged onto the yield surface, but at the wrong point.

Found while running the Sutcliffe et al. (2026) §4.3 von Mises collapse with the DR solver.
Measured with a scratch test (trial deviator with shear, perfect plasticity): final/trial
deviator ratios were −0.65 (normal) vs +0.34 (shear) in 2D and −0.74 vs +0.27 in 3D. The normal
deviator even flipped sign. After the fix, all components scale by the same factor (radial),
and the return converges in one correction, as it must for perfect plasticity.

Effect on the DR collapse (2² MPs/cell, 10 load steps): total DR iterations 188k → 56k, run-out
unchanged (2.2 m). The explicit path's von Mises results change too. A VM `slump_problem` run
completes without NaNs. Drucker–Prager (`DP.jl`) has its own return mapping and was not affected.
