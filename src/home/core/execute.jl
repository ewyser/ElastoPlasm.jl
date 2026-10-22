export execute, execute!

# the plot saved into `folder` after a workflow, when `solver.plot.status`
function save_workflow_plot(p::Problem, workflow!, folder)
    p.solver.plot.status || return nothing
    solver   = p.solver
    quantity = join([plot_entry(k).name for k ∈ solver.plot.what], "_")
    name     = "$(solver.solution)_$(length(p.mesh.prprt.L))d_$(solver.basis.which)_$(solver.basis.trsfr)_" *
               "$(solver.material.elastic)_$(workflow!)_$(quantity).png"
    save_plot((; file=joinpath(mkpath(folder), replace(name, " " => "_"))))
end

# a simulation's HDF5 output, `<folder>/<name>_output.h5`, and the text of its case TOML, if one was written
output_file(sim::Simulation) = joinpath(sim.folder, "$(basename(sim.folder))_output.h5")
case_text(sim::Simulation) = (f = joinpath(sim.folder, "$(basename(sim.folder))_case.toml"); isfile(f) ? read(f, String) : nothing)

"""
    execute!(sim::Simulation; workflows=case_workflows(sim), output=output_file(sim), append=false) -> Simulation
    execute!(file; workflows=case_workflows(sim)) -> file
    execute(sim_or_file; workflows, output, append) -> Simulation

Run a simulation's workflows in order, each on its problem (`workflow!(problem)`): its own
(`sim.workflow`, in `sim.lstps` load steps) unless `workflows` are given. Built-in workflows:
`elastodynamic!`, `elastoplastic!`, `elastoquasistatic!`, `thermodynamic!`. When
`solver.plot.status` is on, a plot is saved into the simulation's folder after each workflow.

On a `.jld2` file (`setup_simulation`), the simulation is loaded, run, and written back after each
workflow; `load_simulation(file)` then holds the results. `execute` runs on a copy (of `sim`, or of
the file's simulation), leaves its argument untouched and returns the run simulation (its plots and
output are still written into the folder).

With `solver.output.status`, the run is also recorded as an HDF5 time series (`open_output`, `record!`)
into `output`: started afresh, or with `append` continued; `execute!(file)` appends, since it continues
the saved state, while `execute` and an in-memory `execute!` start a new one. Read it with
`output_series`, `output_frame`, `load_output`.
"""
function execute!(sim::Simulation; workflows=case_workflows(sim), output=output_file(sim), append::Bool=false)
    p   = sim.problem
    rec = p.solver.output.status ? open_output(output, p; append, case_text=case_text(sim)) : nothing
    try
        task_local_storage(:output_recorder, rec) do
            for workflow! ∈ workflows
                @info workflow_log(p.solver; msg="$workflow!")
                task_local_storage(:output_workflow, string(workflow!))
                workflow!(p)
                record!(p.mpts, p.mesh, p.time, p.solver)             # the state the workflow ends on
                save_workflow_plot(p, workflow!, sim.folder)
            end
        end
    finally
        isnothing(rec) || close(rec)
    end
    return sim
end
function execute!(file::AbstractString; workflows=nothing)
    sim = load_simulation(file)
    for workflow! ∈ something(workflows, case_workflows(sim))
        execute!(sim; workflows=(workflow!,), append=true)
        save_simulation(sim, file)                           # the state after each workflow
    end
    return file
end
execute(sim::Simulation; kwargs...) = execute!(deepcopy(sim); kwargs...)
execute(file::AbstractString; kwargs...) = execute!(load_simulation(file); kwargs...)
