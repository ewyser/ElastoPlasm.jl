# `dynamic_relaxation`'s density/porosity updates break solid mass conservation

**Status: open** (found while storing the solid mass `mpts.s.m`; dynamic relaxation was
deliberately left out of that change).

The explicit path keeps `ρ·Ω = (1−n₀)ρ₀Ω₀` exactly (`deform.jl`: `Ω = J·Ω₀`,
`ρ = (1−n)ρ₀`, `n = 1−(1−n₀)/J`), so it now reads the constant `mpts.s.m` instead of
recomputing it. `dynamic_relaxation/update.jl` does not preserve that invariant:

- **Finite-strain `update`** (`update.jl:44-45`): `ρ = ρ₀/J`, ignoring porosity. After
  the first update `ρ·Ω = ρ₀Ω₀`, not `(1−n₀)ρ₀Ω₀`: the mass jumps by a factor `1/(1−n₀)`
  (×1.11 at the default `n₀ = 0.1`).
- **Infinitesimal-strain `update`** (`update.jl:14-15`): `ρ = ρ/ΔJ` is fine on its own, but
  `n = 1 − (1−n_prev)/J` divides the *previous* porosity by the *total* `J`, compounding every
  step. That's the same bug the removed `deform_fast` had.

The 4 dynamic-relaxation mass sites (`assembly.jl:5`, `fint.jl:13,55,133`) still compute
`ρ[p]*Ω[p]`, so they see these drifting values. The likely fix is to read `mpts.s.m[p]`, as the
explicit path does, and to make both updates use `ρ = (1−n)ρ₀` with `n = 1−(1−n₀)/J`, like
`deform.jl`. That changes DR results, so it wasn't done alongside the explicit-path change.
