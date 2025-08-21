export thermal,ic_thermal 

function ic_thermal(; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    L,nel = [64.1584,64.1584/4.0],[40,10]
    
    @info "Setting up mesh & material point system for $(length(L))d thermal problem"
    # init & kwargs
    instr = kwargser(:instr,kwargs;dim=length(L))
    paths = set_paths(fid,info.sys.out;interactive=false)  
    # mesh & mpts initial conditions
    mesh  = setup_mesh(instr     ; geom = get_geom(nel,L,instr)     )
    cmpr  = setup_cmpr(mesh                                         )                       
    mpts  = setup_mpts(mesh,cmpr ; geom = get_thermal(mesh,cmpr,instr))
    # time parameters
    time  = setup_time(instr     ; te=100.0,tg=10.0,tep=0.0) 
    # plot initial cohesion field
    ms = 0.4*instr[:plot][:dims][1]/mesh.prprt.nel[1]
    opts = (;
        dims    = instr[:plot][:dims],
        what    = ["T"],
        backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
        tit     = L" t = "*string(round(0.0,digits=1))*" [s]",
        file    = joinpath(paths[:plot],"$(mesh.prprt.dim)d_T.png"),
    )
    get_plot_field(mpts,mesh,opts);save_plot(opts)
    # display summary
    @info ic_log(mesh,mpts,time,instr)
    misc = (;
        prefix = "$(mesh.prprt.dim)d_$(instr[:fwrk][:trsfr])"
    )
    return (;mesh,mpts,cmpr,time),(;instr,paths,misc)
end

function thermal(ic::NamedTuple,cfg::NamedTuple; workflow::String="thermodynamic")
    @info "Explicit solution to thermal problem"; config_plot()
    # forward-euler explicit workflow
    out = elastoplasm(deepcopy(ic), deepcopy(cfg); problem = workflow)
    # return output with success flag
    return out = (; out..., success=true,)
end