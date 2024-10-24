function init_mapsto(dim::Number,trsfr::String) 
    if trsfr == "mUSL"
        if dim == 2
            kernel1 = flip2Dp2n(CPU())
            kernel3 = flip23Dn2p(CPU())
        elseif dim == 3
            kernel1 = flip2Dp2n(CPU())
            kernel3 = flip23Dn2p(CPU())
        end
        kernel2 = forwardEuler(CPU())
        kernel3a = kernel_momentum(CPU())
        kernel3b = kernel_velocity(CPU())
        kernel3c = kernel_displacement(CPU())
        return Dict(:map!  => (;p2n! = kernel1 ,solve! = kernel2 ,n2p! = kernel3, ),
                    :augm! => (;p2n! = kernel3a,solve! = kernel3b,Δu!  = kernel3c,),)
    elseif trsfr == "tpicUSL"
        if dim == 2
            kernel1 = tpic2Dp2n(CPU())
            kernel2 = pic23Dn2p(CPU())
        elseif dim == 3
            kernel1 = tpic3Dp2n(CPU())
            kernel2 = pic23Dn2p(CPU())
        end
        return (;p2n! = kernel1, n2p! = kernel2,)
    else
        return throw(ArgumentError("$(trsfr) is not a supported|valid mapping"))
    end    
end
function mapsto!(mpD,meD,g,Δt,instr,whereto) 
    if whereto == "p-->n"
        # initialize nodal quantities
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping to mesh
        instr[:cairn][:mapsto!][:map!].p2n!(mpD,meD,g; ndrange=mpD.nmp);sync(CPU())
    elseif whereto == "solve"
        # viscous damping
        η      = 0.1
        # initialize
        meD.fn.= 0.0
        meD.an.= 0.0
        meD.vn.= 0.0
        # solve momentum equation on the mesh using backend-agnostic kernel
        instr[:cairn][:mapsto!][:map!].solve!(meD,Δt,η; ndrange=meD.nno[end]);sync(CPU())
    elseif whereto == "p<--n"
        instr[:cairn][:mapsto!][:map!].n2p!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
        if instr[:trsfr] == "mUSL"
            # initialize for DM
            meD.pn.= 0.0
            meD.vn.= 0.0
            # accumulate material point contributions
            instr[:cairn][:mapsto!][:augm!].p2n!(mpD,meD; ndrange=mpD.nmp);sync(CPU())
            # solve for nodal incremental displacement
            instr[:cairn][:mapsto!][:augm!].solve!(meD; ndrange=meD.nno[end]);sync(CPU())
            # update material point's displacement
            instr[:cairn][:mapsto!][:augm!].Δu!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
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