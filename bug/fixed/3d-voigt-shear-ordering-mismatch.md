# `_kirchoff_stress`'s 3D Voigt shear ordering disagreed with every consumer

**Status: fixed**, as a side effect of the tensor port.

The old 3D stress builder used slots 4/5/6 = xy/xz/yz, but every `p2n!`
kernel and `Del`-facing code uses 4/5/6 = yz/xz/xy. `get_voigt` now uses the
consumer convention throughout. Note this means **3D finite-strain results
change** (2D unaffected — one shear slot only); not verified against a
reference solution, only made internally self-consistent. 3D remains
generally undertested (see
`bug/known/test-workflow-smpm-gimpm-grid-crossing-instabilities.md`'s "3D
conformity check" section).
