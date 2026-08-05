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
function elastoplasm(sim::S; workflow::Vector{F} = [elastodynamic!]) where {S <:String, F <: Function}
    jldopen(sim) do file
        # unpack mesh, mpts, cmpr, instr, paths as aliases
        mesh,mpts,cmpr   = file["ic/mesh"]  , file["ic/mpts"]  , file["ic/cmpr"]
        time             = file["ic/time"]
        instr,paths,misc = file["cfg/instr"], file["cfg/paths"], file["cfg/misc"]
        # action
        for (k,solver!) ∈ enumerate(workflow)
            @info elastoplasm_log(instr; msg = "$solver!")
            solver!(mpts,mesh,cmpr,time,instr)
        end
        sleep(1.0)
        # postprocessing
        if instr.plot.status
            names     = [v.mpts.name for v in instr.plot.what if haskey(v, :mpts)]
            file_path = joinpath(paths[:plot], "$(misc.prefix)_$(join(names, "_")).png")
            opts = (; file = file_path, )
            save_plot(opts)
        end
    end
    # return success message
    exit_log("(✓) Done! exiting...\n")
    return (; simulation=sim, success=true,)
end
function elastoplasm!(sim::S; workflow::Vector{F} = [elastodynamic!]) where {S <:String, F <: Function}
    jldopen(sim,"r+") do file
        # unpack mesh, mpts, cmpr, instr, paths as aliases
        mesh,mpts,cmpr   = file["ic/mesh"]  , file["ic/mpts"]  , file["ic/cmpr"]
        time             = file["ic/time"]
        instr,paths,misc = file["cfg/instr"], file["cfg/paths"], file["cfg/misc"]
        # action
        for (k,solver!) ∈ enumerate(workflow)
            @info elastoplasm_log(instr; msg = "$solver!")
            solver!(mpts,mesh,cmpr,time,instr)
        end
        sleep(1.0)
        # postprocessing
        if instr.plot.status
            names     = [v.mpts.name for v in instr.plot.what if haskey(v, :mpts)]
            figname   = "$(misc.prefix)_$(lowercase(join(names, "_"))).png"
            file_path = joinpath(paths[:plot], replace(figname, " " => "_"))
            opts = (; file = file_path, )
            save_plot(opts)
        end
        # update initial conditions in jld2 file
        delete!(file, "ic")
        file["ic/mesh"] = mesh
        file["ic/mpts"] = mpts
        file["ic/cmpr"] = cmpr
        file["ic/time"] = time
    end
    # return success message
    exit_log("(✓) Done! exiting...\n")
    return (; simulation=sim, success=true,)
end











































