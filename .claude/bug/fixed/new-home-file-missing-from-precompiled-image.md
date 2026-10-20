# A new file under `src/home/` was silently missing from the precompiled image

**Status: fixed** (branch `refactor-init-problem-spec`, commit `b983487`).

Symptom: `UndefVarError: check_case not defined` right after adding `src/home/init/case/checks.jl`,
with no `Failed to include` warning. `Base.include(ElastoPlasm, ".../checks.jl")` loaded the file
cleanly, so the file itself was fine.

Root cause: `superInc`/`collect_and_include_jls` (`src/boot/include.jl`) walk the directory tree at
precompile time. Julia's cache invalidation tracks only the files that were `include`d, so adding a
file to a scanned directory does not invalidate the image: later sessions loaded the cached image
without the new file's definitions. Earlier new files of the same session (`case/points.jl`,
`case/toml.jl`) loaded only because `boot.jl` was edited in the same step, which did invalidate it.
Confirmed: `isdefined(ElastoPlasm, :check_case)` was `false`, and `true` after touching any included
file.

Fix: `superDir` calls `include_dependency(DIR)` on every directory it walks; a directory's mtime
changes when a file is added or removed. Checked both ways: a probe file added under
`src/home/init/case/` was defined at the next `using`, and gone after deleting it.

Related, documented in `CLAUDE.md` ("Project shape"): the `@warn "Failed to include..."` of a
*broken* file appears only while the package precompiles, not in later sessions that load the
cached image; `include` the file directly to see the error again.
