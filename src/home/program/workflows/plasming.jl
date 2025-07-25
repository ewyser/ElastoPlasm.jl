function plasming!(mp,mesh,cmp,g,T,te,tg,instr)
    @info plasming_logs(instr)
    t,Δt,it,ηmax,ηtot = 0.0,instr[:plot][:freq],0,0,0
    checks = sort(unique([collect(t+Δt:Δt:T);te;T]))
    # action
    prog = Progress(length(checks);dt=0.5,desc="Plasming...",barlen=10)
    for (k,time) ∈ enumerate(checks)
        # plot/save
        savlot(mp,mesh,t,instr)
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
        # update progress bas
        next!(prog;showvalues = get_vals(mesh,mp,it,ηmax,ηtot,t/T,"(✗)"))
    end
    finish!(prog)
    return savlot(mp,mesh,t,instr)
end
export plasming!