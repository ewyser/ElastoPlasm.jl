export elastoplastic!,elastodynamic!,thermodynamic!              

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
function elastodynamic!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D,S}) where {T1,T2,D,E,R,S<:ExplicitSolver}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.te)))
        # action
        # action
        # action
    prog = Progress(length(checks);dt=0.5,desc="Solving elastodynamic...",barlen=10)
    for T ∈ checks
        while time.t[1] < T
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            g,dt = get_spacetime(mpts,mesh,cmpr,time,T)
            # mpm cycle
            ignite(mpts,mesh,instr)
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
function elastoplastic!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D,S}) where {T1,T2,D,E,R,S<:ExplicitSolver}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.t[2])))
    g         = get_g(mesh.prprt; G = T2(9.81))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Solving elastoplastic...",barlen=10)
    for T ∈ checks
        while time.t[1] < T
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(mpts,mesh.prprt,cmpr,time,T)
            # mpm cycle
            ignite(mpts,mesh,instr)
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
function thermodynamic!(mpts::Point{T1,T2,D,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,time::Time{T1,T2},instr::Instruction{T1,T2,D,S}) where {T1,T2,D,E,R,S<:ExplicitSolver}
    it,checks = T1(0), T2.(sort(collect(time.t[1]:instr.plot.freq:time.t[2])))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Solving thermodynamic!...",barlen=10)
    for T ∈ checks
        while time.t[1] < T
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(mpts,mesh.prprt,cmpr,time,T)
            # mpm cycle
            ignite(mpts,mesh     ,instr)
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