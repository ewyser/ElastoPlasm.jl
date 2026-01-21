export elastoplasm,elastoplasm!,elastoplastic!,elastodynamic!,thermodynamic!              

"""
    elastodynamic!(mpts::Point{T1,T2,E,R}, mesh, cmpr::NamedTuple, time::NamedTuple, instr::NamedTuple)

Run the explicit elastodynamic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,E,R}`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `time::NamedTuple`: Time stepping configuration.
- `instr::NamedTuple`: Simulation instructions and options.

# Behavior
- Advances the simulation in time using an explicit MPM cycle, updating material points and mesh.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function elastodynamic!(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,E,R,D}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.te)))
        # action
        # action
        # action
    prog = Progress(length(checks);dt=0.5,desc="Solving elastodynamic...",barlen=10)
    for T ∈ checks
        while T > time.t[1]
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            g,dt = get_spacetime(mpts,mesh,cmpr,time,T)
            # mpm cycle
            shpfun(mpts,mesh,instr)
            mapsto(mpts,mesh.s,g,dt,instr)    
            elasto(mpts,mesh,cmpr,dt,instr)
            # update sim parameters
            time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)
        end
        # plot/save
        savlot(mpts,mesh,time.t[1],instr)
        # update progress bar
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end  
"""
    elastoplastic!(mpts::Point{T1,T2,E,R}, mesh, cmpr::NamedTuple, time::NamedTuple, instr::NamedTuple)

Run the explicit elastoplastic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,E,R}`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `time::NamedTuple`: Time stepping configuration.
- `instr::NamedTuple`: Simulation instructions and options.

# Behavior
- Advances the simulation in time using an explicit MPM cycle with elastoplastic update.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function elastoplastic!(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,E,R,D}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.t[2])))
    g         = get_g(mesh.prprt; G = T2(9.81))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Solving elastoplastic...",barlen=10)
    for T ∈ checks
        while T > time.t[1]
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(mpts,mesh.prprt,cmpr,time,T)
            # mpm cycle
            shpfun(mpts,mesh,instr)
            mapsto(mpts,mesh.s,g,dt,instr)    
            elastoplast(mpts,mesh,cmpr,dt,instr)
            # update sim parameters
            time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)
        end
        # plot/save
        savlot(mpts,mesh,time.t[1],instr)
        # update progress bar
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end  
"""
    thermodynamic!(mpts::Point{T1,T2,E,R},mesh::MeshThermalPhase{T1,T2,D},cmpr::NamedTuple,time::NamedTuple,instr::NamedTuple) where {D}

Run the explicit thermodynamic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,E,R}`: Material point data structure.
- `mesh::MeshThermalPhase{T1,T2,D}`: Mesh data structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `time::NamedTuple`: Time stepping configuration.
- `instr::NamedTuple`: Simulation instructions and options.

# Behavior
- Advances the simulation in time using an explicit MPM cycle with thermodynamic update.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function thermodynamic!(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,E,R,D}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.t[2])))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Solving thermodynamic!...",barlen=10)
    for T ∈ checks
        while T > time.t[1]
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(mpts,mesh.prprt,cmpr,time,T)
            # mpm cycle
            shpfun(mpts,mesh     ,instr)
            mapsto(mpts,mesh.t,dt,instr)    
            thermo(mpts,mesh.t,instr)
            # update sim parameters
            time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)
        end
        # plot/save
        savlot(mpts,mesh,time.t[1],instr)
        # update progress bar
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end  


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
            opts = (;
                file = joinpath(paths[:plot],"$(misc.prefix)_$(join(last.(instr.plot.what),"_")).png"),
            );save_plot(opts)
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
        if instr.plot:status
            opts = (;
                file = joinpath(paths[:plot],"$(misc.prefix)_$(join(last.(instr.plot.what),"_")).png"),
            );save_plot(opts)
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











































