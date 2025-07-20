function init_mapsto(dim::Number,trsfr::String) 
    kernel2 = euler(CPU())
    if trsfr == "musl"
        if dim == 2
            kernel1 = flip_2d_p2n(CPU())
            kernel3 = flip_nd_n2p(CPU())
        elseif dim == 3
            kernel1 = flip_3d_p2n(CPU())
            kernel3 = flip_nd_n2p(CPU())
        end
        kernel3a = augm_momentum(CPU())
        kernel3b = augm_velocity(CPU())
        kernel3c = augm_displacement(CPU())
        return Dict(:map  => (;p2n! = kernel1 ,solve! = kernel2 ,n2p! = kernel3, ),
                    :augm => (;p2n! = kernel3a,solve! = kernel3b,Δu!  = kernel3c,),)
    elseif trsfr == "tpic"
        if dim == 2
            kernel1 = tpic_2d_p2n(CPU())
            kernel3 = pic_nd_n2p(CPU())
        elseif dim == 3
            kernel1 = tpic_3d_p2n(CPU())
            kernel3 = pic_nd_n2p(CPU())
        end
        return Dict(:map  => (;p2n! = kernel1 ,solve! = kernel2 ,n2p! = kernel3, ),)
    else
        return throw(ArgumentError("$(trsfr) is not a supported|valid mapping scheme"))
    end    
end
function mapsto(mpD,meD,g,Δt,instr) 
    # maps material point to node
        p2n(mpD,meD,g,Δt,instr)
    # solve Eulerian momentum equation
        solve(meD,Δt,instr)
    # maps back solution to material point
        n2p(mpD,meD,Δt,instr)
        if instr[:fwrk][:trsfr] == "musl"
            augm(mpD,meD,Δt,instr)
        end
    return nothing
end



















#=
@views function mapstoN!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.oobf.= 0.0
    # mapping back to mesh
    for dim ∈ 1:meD.nD
        lk = ReentrantLock()
        @threads for p ∈ 1:mpD.nmp
            # accumulation
            lock(lk) do 
                if dim == 1 
                    meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p] 
                end
                meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*mpD.v[p,dim])
                meD.oobf[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*g[dim]      )
                meD.oobf[mpD.p2n[:,p],dim].-= mpD.V[p].*(mpD.B[dim:meD.nD:end,:,p]*mpD.σ[:,p]) 
            end
        end
    end
    return nothing
end
=#