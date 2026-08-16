export elastoplasm,elastoplasm!      

"""
    elastoplasm(ic::NamedTuple, cfg::NamedTuple; mode::String="elastodynamic") -> NamedTuple

Run the main simulation workflow for the given initial conditions and configuration.

# Arguments
- `ic::NamedTuple`: Initial conditions (mesh, mpts, cmpr, time).
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
        # unpack mesh, mpts, basis, cmpr, instr, paths as aliases
        mesh,mpts,basis,cmpr = file["ic/mesh"]   , file["ic/mpts"]  , file["ic/basis"], file["ic/cmpr"]
        time                 = file["ic/time"]
        solver,paths,misc    = file["cfg/solver"], file["cfg/paths"], file["cfg/misc"]
        # action
        for (k,workflow!) ∈ enumerate(workflows)
            @info elastoplasm_log(solver; msg = "$workflow!")
            workflow!(mpts,mesh,basis,cmpr,time,solver)
            # postprocessing
            if solver.plot.status
                dimension   = string(Base.unwrap_unionall(typeof(solver)).parameters[3])
                solution    = string(nameof(typeof(solver)))
                deformation = string(solver.fwrk.deform)
                workflow    = string(workflow!)
                quantity    = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
                name        = "$(dimension)_$(solution)_$(deformation)_$(workflow)_$(quantity).png"
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
        # unpack mesh, mpts, basis, cmpr, instr, paths as aliases
        mesh,mpts,basis,cmpr = file["ic/mesh"]   , file["ic/mpts"]  , file["ic/basis"], file["ic/cmpr"]
        time                 = file["ic/time"]
        solver,paths,misc    = file["cfg/solver"], file["cfg/paths"], file["cfg/misc"]
        # action
        for (k,workflow!) ∈ enumerate(workflows)
            @info elastoplasm_log(solver; msg = "$workflow!")
            workflow!(mpts,mesh,basis,cmpr,time,solver)
            # postprocessing
            if solver.plot.status
                dimension   = string(Base.unwrap_unionall(typeof(solver)).parameters[3])
                solution    = string(nameof(typeof(solver)))
                deformation = string(solver.fwrk.deform)
                workflow    = string(workflow!)
                quantity    = join([v.mpts.name for v in solver.plot.what if haskey(v, :mpts)], "_")
                name        = "$(dimension)_$(solution)_$(deformation)_$(workflow)_$(quantity).png"
                path        = joinpath(paths[:plot],replace(name, " " => "_"))
                opts = (; file = path, )
                save_plot(opts)
            end
            # update initial conditions in jld2 file
            delete!(file, "ic")
            file["ic/mesh"]  = mesh
            file["ic/mpts"]  = mpts
            file["ic/basis"] = basis
            file["ic/cmpr"]  = cmpr
            file["ic/time"]  = time
        end
        sleep(1.0)
    end
    # return success message
    exit_log("(✓) Done! exiting...\n")
    return (; simulation=sim, success=true,)
end
















