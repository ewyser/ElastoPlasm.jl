export elastoplasm,elastoplasm!      

"""
    elastoplasm(sim::String; workflows=[elastodynamic!]) -> (; simulation, success)

Run each `workflow!(mpts, mesh, basis, time, solver)` in `workflows`, in order, on the problem
saved at `sim` — the `.jld2` path returned by a setup function such as `slump_problem`.
Built-in workflows: `elastodynamic!`, `elastoplastic!`, `elastoquasistatic!`,
`thermodynamic!`. When `solver.plot.status` is on, a plot is saved after each workflow.

Opens `sim` read-only: the post-run state is discarded. Use `elastoplasm!` to keep it.
"""
function elastoplasm(sim::S; workflows::Vector{F} = [elastodynamic!]) where {S <:String, F <: Function}
    jldopen(sim) do file
        # unpack mesh, mpts, basis, instr, paths as aliases
        problem              = file["ic/problem"]
        mesh,mpts,time       = problem.mesh, problem.mpts, problem.time
        basis                = file["ic/basis"]
        solver,paths,misc    = file["cfg/solver"], file["cfg/paths"], file["cfg/misc"]
        # action
        for (k,workflow!) ∈ enumerate(workflows)
            @info elastoplasm_log(solver; msg = "$workflow!")
            workflow!(mpts,mesh,basis,time,solver)
            # postprocessing
            if solver.plot.status
                dimension   = string(Base.unwrap_unionall(typeof(solver)).parameters[3])
                basisname   = solver.basis.which
                transfer    = solver.basis.trsfr
                solution    = solver.solution
                deformation = solver.material.elastic
                workflow    = string(workflow!)
                quantity    = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
                name        = "$(solution)_$(dimension)d_$(basisname)_$(transfer)_$(deformation)_$(workflow)_$(quantity).png"
                path        = joinpath(paths[:plot],replace(name, " " => "_"))
                opts = (; file = path, )
                save_plot(opts)
            end
        end
        sleep(1.0)
    end
    # return success message
    exit_log("(✓) Done! exiting...\n")
    return (; simulation=sim, success=true,)
end
"""
    elastoplasm!(sim::String; workflows=[elastodynamic!]) -> (; simulation, success)

Same as `elastoplasm`, but writes the post-run `ic/problem` and `ic/basis` back into `sim`
after each workflow, so the results can be loaded and inspected afterwards.
"""
function elastoplasm!(sim::S; workflows::Vector{F} = [elastodynamic!]) where {S <:String, F <: Function}
    jldopen(sim,"r+") do file
        # unpack mesh, mpts, basis, instr, paths as aliases
        problem              = file["ic/problem"]
        mesh,mpts,time       = problem.mesh, problem.mpts, problem.time
        basis                = file["ic/basis"]
        solver,paths,misc    = file["cfg/solver"], file["cfg/paths"], file["cfg/misc"]
        # action
        for (k,workflow!) ∈ enumerate(workflows)
            @info elastoplasm_log(solver; msg = "$workflow!")
            workflow!(mpts,mesh,basis,time,solver)
            # postprocessing
            if solver.plot.status
                dimension   = string(Base.unwrap_unionall(typeof(solver)).parameters[3])
                basisname   = solver.basis.which
                transfer    = solver.basis.trsfr
                solution    = solver.solution
                deformation = solver.material.elastic
                workflow    = string(workflow!)
                quantity    = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
                name        = "$(solution)_$(dimension)d_$(basisname)_$(transfer)_$(deformation)_$(workflow)_$(quantity).png"
                path        = joinpath(paths[:plot],replace(name, " " => "_"))
                opts = (; file = path, )
                save_plot(opts)
            end
            # update initial conditions in jld2 file
            delete!(file, "ic")
            file["ic/problem"] = MechanicalProblem(mesh,mpts,time)
            file["ic/basis"]   = basis
        end
        sleep(1.0)
    end
    # return success message
    exit_log("(✓) Done! exiting...\n")
    return (; simulation=sim, success=true,)
end
















