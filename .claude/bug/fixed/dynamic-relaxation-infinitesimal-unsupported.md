# `dynamic_relaxation`'s plain (non-u-P) path only actually supports `strain.deform="finite"`

**Status: fixed** on `feat-implicit-dr-sutcliffe`.

Despite a comment in `implicit.jl` claiming the opposite: `pt_solve!`/`relax`
always call `oobf_assembly!` with the 4-trailing-arg `ST<:LogarithmicStrain`
(finite-strain) shape; under `"infinitesimal"` only a 3-arg method exists,
throwing `MethodError`. Confirmed pre-existing, not caused by any recent
refactor. `test_column.jl` and every other verified `dynamic_relaxation` run
use the default `"finite"`, which is why this went unnoticed. Someone needs
to either fix `oobf_assembly`'s infinitesimal method to actually get called,
or make `"infinitesimal"` genuinely unsupported by construction (an explicit
error) rather than an opaque `MethodError` deep in a kernel.

## Resolution

Made unsupported by construction: `elastoquasistatic!` throws an `ArgumentError` unless the
strain type is `LogarithmicStrain`. The Sutcliffe et al. (2026) DR scheme it now implements is
finite-strain only, and the old infinitesimal/1D/3D `oobf_assembly` methods were removed.
