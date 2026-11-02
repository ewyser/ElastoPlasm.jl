export elastoplastic!,elastodynamic!,thermodynamic!              

"""
    elastodynamic!(mpts::Point{T1,T2,CM}, mesh, time::NamedTuple, instr::NamedTuple)

Run the explicit elastodynamic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,CM}`: Material point data structure.
- `mesh`: Mesh data structure.
- `time::NamedTuple`: Time stepping configuration.
- `instr::NamedTuple`: Simulation instructions and options.

# Behavior
- Advances the simulation in time using an explicit MPM cycle, updating material points and mesh.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function elastodynamic!(mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},time::Time{T1,T2},solver::ExplicitSolver{T1,T2}) where {T1,T2}
    it,stops = T1(0), checkpoints(time,solver,time.te)     # (t, plot, record): get.jl
    # the arrays on the solver's device (the same objects on the host: device.jl)
    dm,dg,db = to_device((mpts,mesh,basis),solver)
    # action
    prog = Progress(length(stops);dt=0.5,desc="Solving elastodynamic...",barlen=10)
    for (T, plotted, recorded) ∈ stops
        while time.t[1] < T
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            g,dt = get_spacetime(dm,dg,time,T)
            # mpm cycle
            ignite(dm,dg,db,solver)
            mapsto(dm,dg,db,g,dt,solver)
            elasto(dm,dg,db,dt,solver)
            # update sim parameters
            time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)
        end
        # plot, record (from the host problem, brought up to date)
        (plotted || recorded) && to_host!((mpts,mesh,basis),(dm,dg,db))
        plotted  && bake(mpts,mesh,time.t[1],solver)
        recorded && record!(mpts,mesh,time,solver; step=it)
        # update progress bar
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    to_host!((mpts,mesh,basis),(dm,dg,db))
    finish!(prog)
    return nothing
end  
"""
    elastoplastic!(mpts::Point{T1,T2,CM}, mesh, time::NamedTuple, instr::NamedTuple)

Run the explicit elastoplastic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,CM}`: Material point data structure.
- `mesh`: Mesh data structure.
- `time::NamedTuple`: Time stepping configuration.
- `instr::NamedTuple`: Simulation instructions and options.

# Behavior
- Advances the simulation in time using an explicit MPM cycle with elastoplastic update.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function elastoplastic!(mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},time::Time{T1,T2},solver::ExplicitSolver{T1,T2}) where {T1,T2}
    it,stops = T1(0), checkpoints(time,solver,time.t[2])     # (t, plot, record): get.jl
    g         = get_g(mesh.prprt; G = T2(9.81))
    # the arrays on the solver's device (the same objects on the host: device.jl)
    dm,dg,db  = to_device((mpts,mesh,basis),solver)
    # action
    prog = Progress(length(stops);dt=0.5,desc="Solving elastoplastic...",barlen=10)
    for (T, plotted, recorded) ∈ stops
        while time.t[1] < T
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(dm,dg.prprt,time,T)
            # mpm cycle
            ignite(dm,dg,db,solver)
            mapsto(dm,dg,db,g,dt,solver)
            elastoplast(dm,dg,db,dt,solver)
            # update sim parameters
            time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)
        end
        # plot, record (from the host problem, brought up to date)
        (plotted || recorded) && to_host!((mpts,mesh,basis),(dm,dg,db))
        plotted  && bake(mpts,mesh,time.t[1],solver)
        recorded && record!(mpts,mesh,time,solver; step=it)
        # update progress bar
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    to_host!((mpts,mesh,basis),(dm,dg,db))
    finish!(prog)
    return nothing
end  
"""
    thermodynamic!(mpts::Point{T1,T2,D,CM},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},time::Time{T1,T2},instr::S) where {T1,T2,D,CM,S<:Solver{T1,T2,D}}

Run the explicit thermodynamic workflow for the given mesh, material points, constitutive model, and simulation configuration.

# Arguments
- `mpts::Point{T1,T2,D,CM}`: Material point data structure.
- `mesh::Mesh{T1,T2,D}`: Mesh data structure.
- `basis::Basis{T1,T2,D}`: Topology container linking `mesh` and `mpts`.
- `time::Time{T1,T2}`: Time stepping configuration.
- `instr::S`: Solver instance (e.g. `ExplicitSolver`) with matching `T1,T2,D`.

# Behavior
- Advances the simulation in time using an explicit MPM cycle with thermodynamic update.
- Plots and saves results at specified intervals.
- Displays a progress bar.

# Returns
- `nothing`
"""
function thermodynamic!(mpts::Point{T1,T2,D,CM},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},time::Time{T1,T2},instr::S) where {T1,T2,D,CM,S<:Solver{T1,T2,D}}
    host_only(instr, "thermodynamic!")
    it,stops = T1(0), checkpoints(time,instr,time.t[2])     # (t, plot, record): get.jl
    # action
    prog = Progress(length(stops);dt=0.5,desc="Solving thermodynamic!...",barlen=10)
    for (T, plotted, recorded) ∈ stops
        while time.t[1] < T
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(mpts,mesh.prprt,time,T)
            # mpm cycle
            ignite(mpts,mesh,basis,instr)
            mapsto(mpts,mesh.t,basis,dt,instr)
            thermo(mpts,mesh.t,basis,instr)
            # update sim parameters
            time.t[1],it,toc = time.t[1]+dt,it+T1(1),(time_ns()-tic)
        end
        # plot, record
        plotted  && bake(mpts,mesh,time.t[1],instr)
        recorded && record!(mpts,mesh,time,instr; step=it)
        # update progress bar
        next!(prog;showvalues = get_vals(mesh,mpts,it))
    end
    finish!(prog)
    return nothing
end

# the workflows on a whole problem (`execute!(sim)`)
elastodynamic!(p::Problem) = elastodynamic!(p.mpts, p.mesh, p.basis, p.time, p.solver)
elastoplastic!(p::Problem) = elastoplastic!(p.mpts, p.mesh, p.basis, p.time, p.solver)
thermodynamic!(p::Problem) = thermodynamic!(p.mpts, p.mesh, p.basis, p.time, p.solver)
