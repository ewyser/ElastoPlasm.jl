# `MeshSolidPhase.Mᵢⱼ::Matrix{T2}` allocated a full dense `nno×nno` matrix unconditionally in `setup_mesh.jl`, regardless of solver

**Status: fixed (removed).**

O(nodes²) memory: harmless at small mesh sizes (e.g. 6,601 nodes → ~349 MB)
but `OutOfMemoryError` at moderate resolution (32,841 nodes → ~8.6 GB) —
found while trying to run `collapse_problem` at a MaterialPointSolver.jl-
matched fine resolution (`nel=[368,88]`). Confirmed dead (zero reads
anywhere in `src/`, only its own declaration/construction) before removing —
the explicit solver path already uses the separate lumped/diagonal
`m::Vector{T2}` nodal mass field for everything; `Mᵢⱼ` looks like leftover
scaffolding for a never-finished consistent-mass-matrix formulation. Removed
the field from `MeshSolidPhase` (`src/boot/needs/types/problem/eulerian.jl`)
and its construction in `build_solid_mesh_phase` (`src/home/init/setup_mesh.jl`).
