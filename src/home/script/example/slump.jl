export slump,ic_slump 

function ic_slump(L::Vector{Float64},nel::Vector{Int64}; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    @info "Setting up mesh & material point system for $(length(L))d slump problem"
    # init & kwargs
    instr   = kwargser(:instr,kwargs;dim=length(L))
    paths   = set_paths(fid,info.sys.out;interactive=false)      
    # mesh & mp initial conditions
    ni      = 2  
    mesh    = setup_mesh(nel,L,instr)    
    cmpr    = setup_cmpr(mesh.dim,instr)                       
    mp      = setup_mps(mesh,cmpr;define=inislump(mesh,cmpr,ni,instr))
    return (;mesh,mp,cmpr),(;instr,paths)
end

function slump(ic::NamedTuple,cfg::NamedTuple;)
    @info "Explicit solution to slump problem";configPlot()
    # extract mesh, mp, cmpr, instr, paths
    mesh,mp,cmpr = deepcopy(ic[:mesh]  ), deepcopy(ic[:mp]    ), deepcopy(ic[:cmpr])
    instr,paths  = deepcopy(cfg[:instr]), deepcopy(cfg[:paths])
    # plot initial cohesion field
    plotcoh(mp,cmpr,paths)
    # independant physical constant
    g,tg = 9.81,15.0/1.5                                                   
    # constitutive model & time parameters
    T,te = 15.0,10.0                           
    # action
    out  = plasming!(mp,mesh,cmpr,g,T,te,tg,instr)
    # postprocessing
    @info "Fig(s) saved at $(paths[:plot])"
    path =joinpath(paths[:plot],"$(mesh.dim)d_$(mp.nmp)_$(mesh.nel[end])_$(join(instr[:plot][:what]))_$(instr[:basis][:which])_$(instr[:fwrk][:deform])_$(instr[:fwrk][:trsfr])_$(instr[:fwrk][:locking])_$(cmpr[:cmType])_$(instr[:perf])_$(first(instr[:nonloc])).png")
    savefig(path)
    msg("(✓) Done! exiting...\n")
    return sucess=true
end