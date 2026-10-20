export save_problem, load_problem

"""
    save_problem(file, problem; paths) -> file
    load_problem(file) -> Problem

Write a `Problem` to a `.jld2` file, and read it back. The file holds `"problem"` (the `Problem`),
`"paths"` (the output folders, `mkpaths`) and `"misc"`; `elastoplasm`/`elastoplasm!` run what it holds.
"""
function save_problem(file::AbstractString, problem::Problem; paths)
    @info ic_log(problem.mesh, problem.mpts, problem.time, problem.solver)
    misc = (; prefix="$(length(problem.mesh.prprt.nel)-1)d_$(problem.solver.basis.trsfr)")
    jldopen(file, "w") do f
        f["problem"] = problem
        f["paths"]   = paths
        f["misc"]    = misc
    end
    return file
end
function load_problem(file::AbstractString)
    return load(file, "problem")
end
