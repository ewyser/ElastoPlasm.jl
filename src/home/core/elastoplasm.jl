export elastoplasm,elastoplasm!      

"""
    elastoplasm(ic::NamedTuple, cfg::NamedTuple; mode::String="elastodynamic") -> NamedTuple

Run the main simulation workflow for the given initial conditions and configuration.

# Arguments
- `ic::NamedTuple`: Initial conditions (mesh, mpts, basis, time).
- `cfg::NamedTuple`: Simulation configuration (instr, paths).
- `mode::String`: (Optional) Workflow mode: "elastodynamic", "elastoplastic", or "all-in-one" (default: "elastodynamic").

# Behavior
- Runs the selected workflow problem, logging progress and saving results.
- Handles postprocessing and output file naming.
- Returns a named tuple with the input initial conditions and configuration.

# Returns
- `NamedTuple`: Contains the input `ic` and `cfg`.
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
                solution    = solver.solution
                deformation = string(solver.strain.deform)
                workflow    = string(workflow!)
                quantity    = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
                name        = "$(solution)_$(dimension)d_$(basisname)_$(deformation)_$(workflow)_$(quantity).png"
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
                solution    = solver.solution
                deformation = string(solver.strain.deform)
                workflow    = string(workflow!)
                quantity    = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
                name        = "$(solution)_$(dimension)d_$(basisname)_$(deformation)_$(workflow)_$(quantity).png"
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
















