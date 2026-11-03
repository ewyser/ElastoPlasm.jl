# CI on macOS: `KeyError: :dev1` in `get_host`, and the banner prints `vnothing`

**Status: fixed** (branch `refactor-unifiedbackend-submodule`): the host registration by the library-wide fix, in
UnifiedBackend.jl (commit `e931360`), which ElastoPlasm now takes its backends from; the banner by
`pkgversion(@__MODULE__)`.

Seen in the PR #184 CI log on `macOS-latest`. Ubuntu and Windows are not affected by the first
issue.

- **`KeyError: key :dev1 not found`** (`src/boot/needs/backend.jl:92`, `get_host`), before any
  simulation runs. `add_backend!` registers a CPU only if the CPU model string contains a listed
  brand: `"Intel(R)"`/`"AMD"` for the x86_64 entry, `"Apple"`/`"AMD"` for aarch64. CI installs
  `arch: x64` Julia on the Apple Silicon runner, so Julia runs under Rosetta with
  `Sys.ARCH == :x86_64`, while the model string names an Apple chip. No brand matches,
  `bckd.cpu` stays empty, and every `get_solver` call fails. Any CPU whose name lists none of the
  brands (e.g. ARM Linux servers) would hit the same failure. Fix options: register the CPU
  whenever the host architecture entry is functional, without the brand check (library-wide
  fix), or set `arch: aarch64` for the macOS job (CI only).
- **`vnothing` in the welcome banner.** `get_version()` (`src/boot/needs/utils.jl:107`) returns
  `Pkg.project().version`, the version of the *active* environment, which under `Pkg.test` is a
  temporary test environment with no version. `pkgversion(@__MODULE__)` (Julia ≥ 1.9) reads the
  package's own `Project.toml`. Cosmetic.

**Resolution.** ElastoPlasm's own copy of the backend code (`src/boot/needs/backend.jl`) is gone; it takes the
execution backends from UnifiedBackend (`lib/UnifiedBackend`, a submodule), whose `add_backend!` now registers the
architecture's hosts whatever the brand, with brand `"unknown"` when the model names none (`host_brand`), the
architecture as name when the model is empty, and `CPU()`/`Array` for an architecture it does not list; there is always
a `:dev1`. Tested by `host_brand` on an Apple, a Grace (`Neoverse-V2`) and an empty model, and by registering an
unlisted architecture (`:riscv64`). Since then the label is a `Platform`'s `vendor` (`host_vendor`, one list of
makers for every architecture, `"Intel"`/`"AMD"`/`"Apple"`), so x86_64 Julia under Rosetta reports `"Apple"`, not
`"unknown"`. ElastoPlasm's `get_version()` (`src/boot/needs/utils.jl`) reads
`pkgversion(@__MODULE__)`, its own version under any active environment.
