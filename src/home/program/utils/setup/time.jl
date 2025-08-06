function setup_time(ϵ::Type=Float64; te=0.0,tg=0.0,tep=0.0)
    time  = (; 
        t = ϵ.([0.0,te+tep]), 
        te = ϵ(te), 
        tg = if tg > te ϵ(te) else ϵ(tg) end, 
        tep = ϵ(tep),
    )
    return time
end