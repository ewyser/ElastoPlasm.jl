#L = [64.1584,12.80]
function slump(L::Vector{Float64},nel::Int64; fid::String=first(splitext(basename(@__FILE__))),kwargs...)
    @info "Execution of slump()"
    configPlot()
    # init & kwargs
    instr   = kwargser(:instr,kwargs;dim=length(L))
    paths   = set_paths(fid,info.sys.out;interactive=false)
    # independant physical constant
    g,ni    = 9.81,2                                             
    # constitutive model & time parameters
    cmp     = setup_cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                 
    # mesh & mp setup
    mesh    = setup_mesh(nel,L,instr)                          
    mp      = setup_mps(mesh,cmp;define=inislump(mesh,cmp,ni,instr))
    # plot initial cohesion field
    plotcoh(mp,cmp,paths)
    # action
    out     = plasming!(mp,mesh,cmp,g,T,te,tg,instr)
    # postprocessing
    sleep(2.5)
    @info "Fig(s) saved at $(paths[:plot])"
    path =joinpath(paths[:plot],"$(length(L))D_$(nel)_$(join(instr[:plot][:what]))_$(instr[:basis][:which])_$(instr[:fwrk][:deform])_$(instr[:fwrk][:trsfr])_$(instr[:fwrk][:locking])_$(cmp[:cmType])_$(instr[:perf])_$(first(instr[:nonloc])).png")
    savefig(path)
    msg("(✓) Done! exiting...\n")
    return instr
end
export slump