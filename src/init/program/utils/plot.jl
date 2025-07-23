@views function plotcoh(mp,cmParam,paths)
    xp,coh,phi = mp.x,mp.c₀,mp.ϕ
    coh0,ϕ0    = cmParam[:c0],cmParam[:ϕ0]


    coh0 = coh0/1e3
    if size(xp,1)==2
        gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
        scatter(xp[1,:],xp[2,:],zcolor=coh./1e3,
        markershape=:circle,
        label="",
        xlabel = L"$x-$direction [m]",
        ylabel = L"$z-$direction [m]",
        title  = L"c_0(x_p)",
        aspect_ratio=1,
        c=:vik,
        clims=(coh0-coh0/2,coh0+coh0/2),
        markerstrokecolor=:auto,
        markerstrokewidth=0,
        ylim=(-10,20),
        )
    elseif size(xp,1)==3
        gr(size=(2.0*250,4*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
        scatter(xp[1,:],xp[2,:],xp[3,:],zcolor=coh./1e3,
        markershape=:circle,
        label="",
        xlabel = L"$x-$direction",
        ylabel = L"$z-$direction",
        title  = L"c_0(x_p)",
        aspect_ratio=:equal,
        c=:vik,
        clims=(coh0-coh0/2,coh0+coh0/2),
        markerstrokecolor=:auto,
        markerstrokewidth=0,
        )
    end
    savefig(joinpath(paths[:plot],"coh0.png"))
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(xp[1,:],xp[2,:],zcolor=phi,
    markershape=:circle,
    label="",
    xlabel = L"$x-$direction",
    ylabel = L"$z-$direction",
    title  = L"\phi_0(x_p)",
    aspect_ratio=:equal,
    c=:vik,
    clims=(ϕ0-ϕ0/5,ϕ0+ϕ0/5),
    markerstrokecolor=:auto,
    markerstrokewidth=0,
    ylim=(-10,20),
    )
    savefig(joinpath(paths[:plot],"phi0.png"))
end