export plasming!
"""
    plasming!(mp, mesh, cmpr, time, instr) -> Any

Runs the explicit time-stepping workflow for ElastoPlasm, updating material points and mesh state over time.

# Arguments
- `mp`: Material point data structure.
- `mesh`: Mesh data structure.
- `cmpr`: Compression or constitutive model data.
- `time`: Time data structure.
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
result = plasming!(mp, mesh, cmpr, g, 10.0, 5.0, 2.0, instr)
```
"""
function plasming!(mp,mesh,cmpr,time,instr)
    @info plasming_logs(instr)
    it,ηmax,ηtot = 0, 0, 0
    # action
    prog = Progress(length(time.checks);dt=0.5,desc="Plasming...",barlen=10)
    for T ∈ time.checks
        while T > time.t[1]
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            g,dt = get_spacetime(mp,mesh,cmpr,instr,time,T)
            # mpm cycle
            shpfun(mp,mesh,instr)
            mapsto(mp,mesh,g,dt,instr)    
            ηmax = elastoplast(mp,mesh,cmpr,dt,instr)
            # update sim parameters
            time.t[1],it   = time.t[1]+dt     ,it+1
            toc      ,ηtot = ((time_ns()-tic)),max(ηmax,ηtot)
        end
        # plot/save
        savlot(mp,mesh,time.t[1],instr)
        # update progress bar
        next!(prog;showvalues = get_vals(mesh,mp,it,ηmax,ηtot))
    end
    finish!(prog); 
    return sleep(1.0)
end                        