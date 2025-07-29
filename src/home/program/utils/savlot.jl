@views function savlot(mp,mesh,t,instr) 
    if instr[:plot][:status]        
        ms = 0.4*instr[:plot][:dims][1]/mesh.nel[1]
        opts = (;
            dims  = instr[:plot][:dims],
            what = instr[:plot][:what],
            backend = gr(legend=true,markersize=ms,markershape=:circle,markerstrokewidth=0.75,),
            tit     = L" t = "*string(round(t,digits=1))*" [s]",
        )
        return get_plot_field(mp,mesh,opts)
    else
        nothing
    end
    return nothing
end