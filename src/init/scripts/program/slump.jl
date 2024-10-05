function slump(L::Vector{Float64},nel::Int64; kwargs...)
    configPlot()
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "init slump geometry"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,cmParam,instr[:GRF],typeD)                      # material point geometry setup
    # action
    out     = ϵp23De!(mpD,meD,cmParam,g,T,te,tg,instr)
    # postprocessing
    sleep(2.5)
    @info "figs saved in" path_plot
    savefig(path_plot*"$(length(L))D_$(last(instr[:plot]))_$(instr[:shpfun])_$(instr[:fwrk])_$(instr[:trsfr])_$(instr[:vollock])_$(cmParam[:cmType]).png")
    return msg("(✓) Done! exiting...")
end
export slump