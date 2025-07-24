@views function get_dt(mp,mesh,yd,t,time)
    if length(mesh.h)==2
        Δx,Δz = mesh.h
        cmax  = [Δx/(mp.vmax[1]+yd) Δz/(mp.vmax[2]+yd)]; mp.vmax.=0.0 
        dt    = 0.5*maximum(cmax)
    elseif length(h)==3
        Δx   = h[1]
        Δy   = h[2]
        Δz   = h[3]
        vmax = [abs.(vp[1,:]) abs.(vp[2,:]) abs.(vp[3,:])]
        cmax = [maximum(vmax[1,:]) maximum(vmax[2,:]) maximum(vmax[3,:])]
        cmax = [Δx/(cmax[1]+yd) Δy/(cmax[2]+yd) Δz/(cmax[3]+yd)]
        dt   = 0.5*maximum(cmax)
    else
        dt = nothing    
    end
    return min(dt,time-t)
end
function get_g(t::Float64,tg::Float64,ndim::Int64)
    g = 0.0
    if t<=tg 
        g = 9.81*t/tg 
    else
        g = 9.81
    end
    return if ndim == 1 g = [-g] elseif ndim == 2 g = [0.0 -g] elseif ndim == 3 g = [0.0 0.0 -g] end
end

function get_spacetime(mp,mesh,cmp,instr,t,tg,te,time) # t = tw 
    # check elastoplastic status
    if t > te 
        instr[:plast][:status] = true 
    end
    # calculte dt
    cmax = mesh.h./(mp.vmax.+cmp[:c]); mp.vmax.=0.0 
    dt   = min(0.5*maximum(cmax),time-t)
    # ramp-up gravity
    if t<=tg 
        g = 9.81*t/tg 
    else
        g = 9.81
    end
    if mesh.dim == 1 
        g = [-g] 
    elseif mesh.dim == 2 
        g = [0.0 -g] 
    elseif mesh.dim == 3 
        g = [0.0 0.0 -g] 
    end
    return g,dt
end