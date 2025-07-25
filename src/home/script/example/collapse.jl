function collapse(dim,nel,ν,E,ρ0,l0; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Execution of collapse()"
    if dim == 2
        L       = [10.0,1.25*l0     ]                                            # domain geometry
    elseif dim == 3
        L       = [10.0,10.0,1.25*l0]                                            # domain geometry
    end
    # independant physical constant
    instr   = kwargser(:instr,kwargs;dim=length(L))
    paths   = set_paths(fid,info.sys.out;interactive=false)
    # independant physical constant
    g,ni    = 9.81,2                                             
    # constitutive model & time parameters
    cmp     = setup_cm(length(L),instr; E=E,ν=ν,ρ0=ρ0)
    tg      = ceil((1.0/cmp.c)*(2.0*l0)*40.0)
    T,te    = 1.25*tg,1.25*tg   
    # mesh & mp setup
    mesh    = setup_mesh(nel,L,instr)    
    setgeom = inicollapse(mesh,cmp,ni;ℓ₀=l0) 
    mp      = setup_mps(mesh,cmp;define=setgeom)                                        
    z0      = copy(mp.x[:,end])
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
export collapse

#=
plot = (;status=true,freq=1.0,what=["P"],dims=(500.0,250.0),)
dim,nel  = 2,[5,5,10]
# initial parameters 
ν,E,ρ0,l0 = 0.0,1.0e4,80.0,10.0

collapse(dim,nel,ν,E,ρ0,l0; plot)
=#
# collapse(2,10,0.3,1.0e6,2700.0,50.0; plot)