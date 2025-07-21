@views function get_Δt(vp,h,yd,tw,t)
    if length(h)==2
        Δx   = h[1]
        Δz   = h[2]
        vmax = [abs.(vp[:,1]) abs.(vp[:,2])]
        cmax = [maximum(vmax[:,1]) maximum(vmax[:,2])]
        cmax = [Δx/(cmax[1]+yd) Δz/(cmax[2]+yd)]
        Δt   = 0.5*maximum(cmax)
    elseif length(h)==3
        Δx   = h[1]
        Δy   = h[2]
        Δz   = h[3]
        vmax = [abs.(vp[:,1]) abs.(vp[:,2]) abs.(vp[:,3])]
        cmax = [maximum(vmax[:,1]) maximum(vmax[:,2]) maximum(vmax[:,3])]
        cmax = [Δx/(cmax[1]+yd) Δy/(cmax[2]+yd) Δz/(cmax[3]+yd)]
        Δt   = 0.5*maximum(cmax)
    else
        Δt = nothing    
    end
    return min(Δt,t-tw)
end
function get_g(tw::Float64,tg::Float64,ndim::Int64)
    g = 0.0
    if tw<=tg 
        g = 9.81*tw/tg 
    else
        g = 9.81
    end
    return if ndim == 1 g = [-g] elseif ndim == 2 g = [0.0 -g] elseif ndim == 3 g = [0.0 0.0 -g] end
end