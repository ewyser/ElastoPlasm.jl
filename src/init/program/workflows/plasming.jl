function plasming!(mp,mesh,cmParam,g,T,te,tg,instr)
    @info """
    launching ϵlastσPlasm 👻 v$(getVersion()):
    - $(nthreads()) active thread(s) 
    - $(instr[:fwrk][:deform]) strain formulation
    - $(instr[:basis][:which]) calculation cycle
    - $(if instr[:fwrk][:locking] "F-bar locking mitigation" else "no locking mitigation" end)
    - $(if instr[:nonloc][:status] "non-local plastic regularization" else nonlocal = "local plastic formulation" end)
    """
    t,tC,it,ηmax,ηtot = 0.0,instr[:plot][:freq],0,0,0
    # action
    
    prog = ProgressUnknown("plasming...", spinner=true,showspeed=true)
    for (k,time) ∈ enumerate(sort(unique([collect(t+tC:tC:T);te;T])))
        if t > te 
            instr[:plast][:status] = true 
        end
        # plot/save
        savlot(mp,mesh,t,instr)
        while t<time
            # set clock on/off
            tic = time_ns()
            # adaptative Δt & linear increase of gravity
            Δt,g  = get_Δt(mp.v,mesh.h,cmParam[:c],t,T),get_g(t,tg,mesh.dim)
            # mpm cycle
            shpfun(mp,mesh,instr)
            mapsto(mp,mesh,g,Δt,instr)    
            ηmax = elastoplast(mp,mesh,cmParam,Δt,instr)
            # update sim time
            t,it,toc,ηtot = t+Δt,it+1,((time_ns()-tic)),max(ηmax,ηtot)
            # update progress bas
            next!(prog;showvalues = getvals(mesh,mp,it,ηmax,ηtot,t/T,"(✗)"))
        end
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = getvals(mesh,mp,it,ηmax,ηtot,1.0,"(✓)"))
    return savlot(mp,mesh,t,instr)
end
export plasming!