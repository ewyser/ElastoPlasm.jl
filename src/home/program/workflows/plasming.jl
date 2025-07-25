export plasming!
"""
    plasming!(mp, mesh, cmp, g, T, te, tg, instr) -> Any

Runs the explicit time-stepping workflow for ElastoPlasm, updating material points and mesh state over time.

# Arguments
- `mp`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmp`: Compression or constitutive model data.
- `g`: Gravity vector or field.
- `T`: Total simulation time.
- `te`: End time for gravity ramp.
- `tg`: Gravity ramp duration.
- `instr`: Named tuple containing simulation instructions and plotting options.

# Returns
- The result of the final call to `savlot`, typically a plot or saved state.

# Behavior
- Initializes simulation time and progress bar.
- Advances the simulation in time steps, updating mesh and material points.
- Handles adaptive time stepping and gravity ramping.
- Plots and saves results at specified intervals.
- Tracks and displays progress.

# Example
```julia
result = plasming!(mp, mesh, cmp, g, 10.0, 5.0, 2.0, instr)
```
"""
function plasming!(mp,mesh,cmp,time,instr)
    @info plasming_logs(instr)
    t,Δt,it,ηmax,ηtot = 0.0,instr[:plot][:freq],0,0,0
    te,tg,T = time.te,time.tg,time.T
    checks  = sort(unique([collect(t:Δt:T);te;T]))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Plasming...",barlen=10)
    for (k,time) ∈ enumerate(checks)
        while t<time
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            g,dt = get_spacetime(mp,mesh,cmp,instr,t,tg,te,time)
            # mpm cycle
            shpfun(mp,mesh,instr)
            mapsto(mp,mesh,g,dt,instr)    
            ηmax = elastoplast(mp,mesh,cmp,dt,instr)
            # update sim time
            t,it,toc,ηtot = t+dt,it+1,((time_ns()-tic)),max(ηmax,ηtot)
        end
        # plot/save
        savlot(mp,mesh,t,instr)
        # update progress bas
        next!(prog;showvalues = get_vals(mesh,mp,it,ηmax,ηtot,t/T,"(✗)"))
    end
    finish!(prog); sleep(1.0)
    return savlot(mp,mesh,t,instr)
end                        
