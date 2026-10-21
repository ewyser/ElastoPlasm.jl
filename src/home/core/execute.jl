export execute, execute!

# the plot saved into `folder` after a workflow, when `solver.plot.status`
function save_workflow_plot(p::Problem, workflow!, folder)
    p.solver.plot.status || return nothing
    solver   = p.solver
    quantity = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
    name     = "$(solver.solution)_$(length(p.mesh.prprt.L))d_$(solver.basis.which)_$(solver.basis.trsfr)_" *
               "$(solver.material.elastic)_$(workflow!)_$(quantity).png"
    save_plot((; file=joinpath(mkpath(folder), replace(name, " " => "_"))))
end

"""
    execute!(sim::Simulation; workflows=case_workflows(sim)) -> Simulation
    execute!(file; workflows=case_workflows(sim)) -> file
    execute(sim_or_file; workflows) -> Simulation

Run a simulation's workflows in order, each on its problem (`workflow!(problem)`): its own
(`sim.workflow`, in `sim.lstps` load steps) unless `workflows` are given. Built-in workflows:
`elastodynamic!`, `elastoplastic!`, `elastoquasistatic!`, `thermodynamic!`. When
`solver.plot.status` is on, a plot is saved into the simulation's folder after each workflow.

On a `.jld2` file (`setup_simulation`), the simulation is loaded, run, and written back after each
workflow; `load_simulation(file)` then holds the results. `execute` runs on a copy (of `sim`, or of
the file's simulation), leaves its argument untouched and returns the run simulation.
"""
function execute!(sim::Simulation; workflows=case_workflows(sim))
    for workflow! ∈ workflows
        @info elastoplasm_log(sim.problem.solver; msg="$workflow!")
        workflow!(sim.problem)
        save_workflow_plot(sim.problem, workflow!, sim.folder)
    end
    return sim
end
function execute!(file::AbstractString; workflows=nothing)
    sim = load_simulation(file)
    for workflow! ∈ something(workflows, case_workflows(sim))
        execute!(sim; workflows=(workflow!,))
        jldsave(file; simulation=sim)                        # the state after each workflow
    end
    return file
end
execute(sim::Simulation; kwargs...) = execute!(deepcopy(sim); kwargs...)
execute(file::AbstractString; kwargs...) = execute!(load_simulation(file); kwargs...)
