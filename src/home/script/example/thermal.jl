export ic_thermal 

function ic_thermal(; fid::String=first(splitext(basename(@__FILE__))), kwargs...)
    L,nel = [16.0,16.0],[80,80]
    
    @info "Setting up mesh & material point system for $(length(L))d thermal problem"
    # init & kwargs
    instr = kwargser(kwargs; dim=length(L))
    instr = merge(instr, 
        (;
            basis = (;
                which = "bsmpm",
                how = nothing,
            ),
        )
    )
    instr = merge(instr, 
        (;
            plot = (;
                status=true,
                freq=0.01,
                dpi=500,
                what=[("mpts","T"),("mesh","T")],
                cblim  = [(0.0, 25.0),(0.0, 25.0),],
            )
        )
    )
    instr = merge(instr, 
        (;
            bcs   = (;
                dirichlet = [
                    :fixed :fixed;
                    :fixed :fixed], # for 2d, this translates to [lower_x upper_x;lower_y upper_y]
            ),
        )
    )
    paths = set_paths(fid,self.sys.out;interactive=false)  
    # mesh & mpts initial conditions
    mesh  = setup_mesh(instr     ; geom = get_geom(nel,L,instr)       )
    mat   = setup_matconst(mesh                                       )
    mpts  = setup_mpts(mesh,instr,mat ; geom = get_thermal(mesh,mat,instr))
    # time parameters
    time  = setup_time(instr     ; te=5.0,tg=5.0,tep=0.0) 
    # plot initial cohesion field
    dims  = instr[:plot][:dpi].*(mesh.prprt.L[1]./mesh.prprt.L)
    ms    = mesh.prprt.nel[1]*dims[1]
    opts = (;
        dims    = instr[:plot][:dpi].*(mesh.prprt.L./mesh.prprt.L[1]), 
        what    = [("mpts","T"),],
        backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
        tit     = L" t = "*string(round(0.0,digits=1))*" [s]",
        cblim   = [(0.0,20.0),],
        xlim    = (minimum(getindex.(mesh.x, 1)), maximum(getindex.(mesh.x, 1))),
        ylim    = (minimum(getindex.(mesh.x, 2)), maximum(getindex.(mesh.x, 2))),
        file    = joinpath(paths[:plot],"$((length(mesh.prprt.nel)-1))d_T.png"),
    )
    get_plot_field(mpts,mesh,opts);save_plot(opts)
    # display summary
    @info ic_log(mesh,mpts,time,instr)
    misc = (;
        prefix = "$((length(mesh.prprt.nel)-1))d_$(instr[:transfer][:trsfr])"
    )
    # export to jld2 file and return path
    return export_setup(mesh,mpts,time,instr,paths,misc; path = paths[:dat], file = "thermal_simulation")
end