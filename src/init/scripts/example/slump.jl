#L = [64.1584,12.80]
function slump(L::Vector{Float64},nel::Int64; kwargs...)
    @info "execution of slump()"
    configPlot()
    # init & kwargs
    instr  = kwargser(:instr,kwargs;dim=length(L))
    fid    = splitext(basename(@__FILE__))
    paths  = setPaths(first(fid), sys.out;interactive=false)
    # independant physical constant
    g       = 9.81   
    ni      = 2                                            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                 
    # mesh & mp setup
    meD     = meshSetup(nel,L,instr)      
    setgeom = inislump(meD,cmParam,ni,instr)                       
    mpD     = pointSetup(meD,cmParam,instr;define=setgeom)
    # plot initial cohesion field
    plotcoh(mpD,cmParam,paths)
    # action
    out     = plasming!(mpD,meD,cmParam,g,T,te,tg,instr)
    # postprocessing
    sleep(2.5)
    @info "fig(s) saved at $(paths[:plot])"
    path =joinpath(paths[:plot],"$(length(L))D_$(nel)_$(join(instr[:plot][:what]))_$(instr[:basis][:which])_$(instr[:fwrk][:deform])_$(instr[:fwrk][:trsfr])_$(instr[:fwrk][:locking])_$(cmParam[:cmType])_$(instr[:perf])_$(first(instr[:nonloc])).png")
    savefig(path)
    msg("(✓) Done! exiting...\n")
    return instr
end
export slump