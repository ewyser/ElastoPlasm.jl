export ϵlastσplasm,ϵlastσplastic!,ϵlastσdynamic!              

function ϵlastσdynamic!(mp,mesh,cmpr,time,instr)
    it,ηmax,ηtot = 0, 0, 0
    checks = sort(collect(time.t[1]:instr[:plot][:freq]:time.te))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Solving elastodynamic...",barlen=10)
    for T ∈ checks
        while T > time.t[1]
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            g,dt = get_spacetime(mp,mesh,cmpr,instr,time,T)
            # mpm cycle
            shpfun(mp,mesh,instr)
            mapsto(mp,mesh,g,dt,instr)    
            elasto(mp,mesh,cmpr,dt,instr)
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
function ϵlastσplastic!(mp,mesh,cmpr,time,instr)
    it,ηmax,ηtot = 0, 0, 0
    checks = sort(collect(time.t[1]:instr[:plot][:freq]:time.t[2]))
    g = get_g(mesh.dim)
    # action
    prog = Progress(length(checks);dt=0.5,desc="Solving elastoplastic...",barlen=10)
    for T ∈ checks
        while T > time.t[1]
            # set clock on/off
            tic = time_ns()
            # adaptative dt & linear increase of gravity
            dt  = get_dt(mp,mesh,cmpr,instr,time,T)
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

function ϵlastσplasm(mp,mesh,cmpr,time,paths,instr; load::String="elastodynamic")
    # action
    if load == "elastodynamic"
        @info elastoplasm_log(instr;           )
        ϵlastσdynamic!(mp,mesh,cmpr,time,instr )
    elseif load == "elastoplastic"
        @info elastoplasm_log(instr; msg = load) 
        ϵlastσplastic!(mp,mesh,cmpr,time,instr )
    else
        @error "Invalid workflow: $(load). Choose 'elastodynamic' or 'elastoplastic'."
        return false
    end
    # postprocessing
    @info "Fig(s) saved at $(paths[:plot])"
    path =joinpath(paths[:plot],"$(mesh.dim)d_$(mp.nmp)_$(mesh.nel[end])_$(join(instr[:plot][:what]))_$(instr[:basis][:which])_$(instr[:fwrk][:deform])_$(instr[:fwrk][:trsfr])_$(instr[:fwrk][:locking])_$(cmpr[:cmType])_$(instr[:perf])_$(first(instr[:nonloc])).png")
    savefig(path)
    # return success message
    return msg("(✓) Done! exiting...\n")
end