function init_mapsto(dim::Number,trsfr::String;what::String="p2n!") 
    if trsfr == "mUSL"
        if what == "p2n!"
            if dim == 2
                p2n! = flip2Dp2n(CPU())
            elseif dim == 3
                p2n! = flip3Dp2n(CPU())
            end
            return p2n!
        elseif what == "n2p!"
            n2p! = flip23Dn2p(CPU())
            return n2p!
        end
    elseif trsfr == "tpicUSL"
        if what == "p2n!"
            if dim == 2
                p2n! = tpic2Dp2n(CPU())
            elseif dim == 3
                p2n! = tpic3Dp2n(CPU())
            end
            return p2n!
        elseif what == "n2p!"
            n2p! = pic23Dn2p(CPU())
            return n2p!
        end
    else
        return throw(ArgumentError("$(trsfr) is not a supported|valid mapping"))
    end    
end
function mapsto!(mpD,meD,g,Δt,instr,whereto) 
    if whereto == "p>n"
        # initialize nodal quantities
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping to mesh
        instr[:cairn].p2n!(mpD,meD,g; ndrange=mpD.nmp);sync(CPU())
    elseif whereto == "p<n"
        instr[:cairn].n2p!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
        if instr[:trsfr] == "mUSL"
            DM!(mpD,meD,Δt,instr)
        end
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