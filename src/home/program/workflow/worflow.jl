export elastoplasm,elastoplastic!,elastodynamic!              

"""
    elastodynamic!(mpts::Point{T1,T2}, mesh, cmpr::NamedTuple, time::NamedTuple, instr::NamedTuple)

Run the explicit elastodynamic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
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
function elastodynamic!(mpts::Point{T1,T2},mesh,cmpr::NamedTuple,time::NamedTuple,instr::NamedTuple) where {T1,T2}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr[:plot][:freq]:time.te)))
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
            mapsto(mpts,mesh,g,dt,instr)    
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
    elastoplastic!(mpts::Point{T1,T2}, mesh, cmpr::NamedTuple, time::NamedTuple, instr::NamedTuple)

Run the explicit elastoplastic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
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
function elastoplastic!(mpts::Point{T1,T2},mesh,cmpr::NamedTuple,time::NamedTuple,instr::NamedTuple) where {T1,T2}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr[:plot][:freq]:time.t[2])))
    g         = get_g(mesh; G = T2(9.81))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Solving elastoplastic...",barlen=10)
    for T ∈ checks
        while T > time.t[1]
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(mpts,mesh,cmpr,time,T)
            # mpm cycle
            shpfun(mpts,mesh,instr)
            mapsto(mpts,mesh,g,dt,instr)    
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
    elastoplasm(ic::NamedTuple, cfg::NamedTuple; mode::String="elastodynamic") -> NamedTuple

Run the main simulation workflow for the given initial conditions and configuration.

# Arguments
- `ic::NamedTuple`: Initial conditions (mesh, mpts, cmpr, time).
- `cfg::NamedTuple`: Simulation configuration (instr, paths).
- `mode::String`: (Optional) Workflow mode: "elastodynamic", "elastoplastic", or "all-in-one" (default: "elastodynamic").

# Behavior
- Runs the selected workflow, logging progress and saving results.
- Handles postprocessing and output file naming.
- Returns a named tuple with the input initial conditions and configuration.

# Returns
- `NamedTuple`: Contains the input `ic` and `cfg`.
"""
function elastoplasm(ic::NamedTuple,cfg::NamedTuple; mode::String="elastodynamic")
    # unpack mesh, mpts, cmpr, instr, paths as aliases
    mesh,mpts,cmpr   = ic[:mesh]  , ic[:mpts]  , (ic[:cmpr])
    time             = ic[:time]
    instr,paths,misc = cfg[:instr], cfg[:paths], cfg[:misc]
    # action
    @info elastoplasm_log(instr; msg = mode) 
    if mode == "elastodynamic"
        elastodynamic!(mpts,mesh,cmpr,time,instr)
    elseif mode == "elastoplastic"
        elastoplastic!(mpts,mesh,cmpr,time,instr)
    elseif mode == "all-in-one"
        elastodynamic!(mpts,mesh,cmpr,time,instr)
        elastoplastic!(mpts,mesh,cmpr,time,instr)
    else
        @error "Invalid workflow: $(mode). Choose 'elastodynamic', 'elastoplastic' or 'all-in-one'."
        return false
    end
    sleep(1.0)
    # postprocessing
    if instr[:plot][:status]
        opts = (;
            file = joinpath(paths[:plot],"$(misc[:file]).png"),
        );save_plot(opts)
    end
    # return success message
    exit_log("(✓) Done! exiting...\n")
    return out = (; ic,cfg,)
end

#file = joinpath(paths[:plot],"$(mesh.dim)d_$(instr[:fwrk][:trsfr])_$(mode).png"),