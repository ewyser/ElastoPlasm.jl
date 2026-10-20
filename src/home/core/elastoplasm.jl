export elastoplasm,elastoplasm!

# a workflow on a problem: `workflow!(problem)` when it takes one (the built-ins, closures such as
# `(a...) -> elastoquasistatic!(a...; lstps=10)`), else `workflow!(mpts, mesh, basis, time, solver)`
run_workflow!(workflow!, p::Problem) = applicable(workflow!, p) ? workflow!(p) : workflow!(p.mpts, p.mesh, p.basis, p.time, p.solver)

# the plot saved after a workflow when `solver.plot.status` (into `paths[:plot]`; none without paths)
function save_workflow_plot(p::Problem, workflow!, paths)
    (p.solver.plot.status && !isnothing(paths)) || return nothing
    solver   = p.solver
    quantity = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
    name     = "$(solver.solution)_$(length(p.mesh.prprt.L))d_$(solver.basis.which)_$(solver.basis.trsfr)_" *
               "$(solver.material.elastic)_$(workflow!)_$(quantity).png"
    save_plot((; file=joinpath(paths[:plot], replace(name, " " => "_"))))
end

# the HDF5 output of a simulation file: `<name>_output.h5` beside `<name>_simulation.jld2`, and its case TOML
output_file(sim::AbstractString) = joinpath(dirname(sim), replace(basename(sim), r"_simulation\.jld2$" => "", r"\.jld2$" => "") * "_output.h5")
case_text(sim::AbstractString) = (f = joinpath(dirname(sim), replace(basename(sim), r"_simulation\.jld2$" => "") * "_case.toml"); isfile(f) ? read(f, String) : nothing)

"""
    elastoplasm(problem::Problem; workflows=[elastodynamic!], paths=nothing, output=nothing, append=false) -> Problem
    elastoplasm(file; workflows=[elastodynamic!]) -> (; simulation, success)
    elastoplasm!(file; workflows=[elastodynamic!]) -> (; simulation, success)

Run each workflow of `workflows`, in order, on a problem: `workflow!(problem)`, or
`workflow!(mpts, mesh, basis, time, solver)` for a workflow written that way. Built-in workflows:
`elastodynamic!`, `elastoplastic!`, `elastoquasistatic!`, `thermodynamic!`. When
`solver.plot.status` is on, a plot is saved after each workflow (into `paths[:plot]`).

`elastoplasm(problem)` runs in memory and returns the (changed) problem. On a `.jld2` file (from
`setup_simulation` or a `*_problem` call), `elastoplasm` runs a copy and discards the result;
`elastoplasm!` writes the problem back after each workflow, so the results can be loaded
(`load_problem`) and inspected afterwards.

With `solver.output.status`, the run is also recorded as an HDF5 time series (`open_output`,
`record!`): into `output` (in memory), or `<name>_output.h5` beside the `.jld2`, which `elastoplasm`
starts afresh and `elastoplasm!` continues (`append`), since it continues from the saved state. Read
it with `output_series`, `output_frame`, `load_output`.
"""
function elastoplasm(p::Problem; workflows::Vector{<:Function}=Function[elastodynamic!], paths=nothing,
                     output=nothing, append::Bool=false, case_text=nothing)
    rec = nothing
    if p.solver.output.status
        isnothing(output) && !isnothing(paths) && (output = joinpath(paths[:dat], "output.h5"))
        isnothing(output) ? @warn("output.status is on but no output file is given (output=…): not recorded") :
                            (rec = open_output(output, p; append, case_text))
    end
    try
        task_local_storage(:elastoplasm_output, rec) do
            for workflow! ∈ workflows
                @info elastoplasm_log(p.solver; msg = "$workflow!")
                task_local_storage(:elastoplasm_workflow, string(workflow!))
                run_workflow!(workflow!, p)
                record!(p.mpts, p.mesh, p.time, p.solver)                # the state the workflow ends on
                save_workflow_plot(p, workflow!, paths)
            end
        end
    finally
        isnothing(rec) || close(rec)
    end
    return p
end
function elastoplasm(sim::String; workflows::Vector{<:Function}=Function[elastodynamic!])
    p, paths = load(sim, "problem", "paths")
    elastoplasm(p; workflows, paths, output=output_file(sim), case_text=case_text(sim))
    sleep(1.0)
    exit_log("(✓) Done! exiting...\n")
    return (; simulation=sim, success=true,)
end
function elastoplasm!(sim::String; workflows::Vector{<:Function}=Function[elastodynamic!])
    p, paths = load(sim, "problem", "paths")
    for workflow! ∈ workflows
        elastoplasm(p; workflows=Function[workflow!], paths, output=output_file(sim), append=true, case_text=case_text(sim))
        jldopen(sim, "r+") do f                               # the post-run state, after each workflow
            delete!(f, "problem")
            f["problem"] = p
        end
    end
    sleep(1.0)
    exit_log("(✓) Done! exiting...\n")
    return (; simulation=sim, success=true,)
end
