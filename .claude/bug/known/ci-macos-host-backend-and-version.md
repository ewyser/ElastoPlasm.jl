# CI on macOS: `KeyError: :dev1` in `get_host`, and the banner prints `vnothing`

**Status: open** (diagnosed; fix deferred by the maintainer).

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
