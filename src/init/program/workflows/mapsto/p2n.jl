@kernel inbounds = true function flip_1d_p2n(mpD,meD,g)
    p = @index(Global)
    if p≤mpD.nmp 
        # accumulation
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            @atom meD.pn[no]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p])
            # lumped mass matrix
            @atom meD.mn[no]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
            # consistent mass matrix
            # meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p]   
            @atom meD.oobf[no]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g) 
            @atom meD.oobf[no]-= mpD.Ω[p]*mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]
        end
    end
end
@kernel inbounds = true function flip_2d_p2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    # lumped mass matrix
                    @atom meD.mn[no]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    # consistent mass matrix
                    # meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p]    
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[3,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p])
                end
            end
        end
    end
end
@kernel inbounds = true function flip_3d_p2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    @atom meD.mn[no      ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p] 
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[5,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[4,p])
                elseif dim == 3
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[3,p])
                end
            end
        end
    end
end
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TPIC transfer scheme, see Nakamura etal, 2023, https://doi.org/10.1016/j.cma.2022.115720
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
@kernel inbounds = true function tpic_2d_p2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]*(mpD.v[p,dim]+mpD.∇vᵢⱼ[dim,1,p]*mpD.δnp[nn,1,p]+mpD.∇vᵢⱼ[dim,2,p]*mpD.δnp[nn,2,p])
                if dim == 1
                    @atom meD.mn[no]      += mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[3,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p])
                end
            end
        end
    end

end
@kernel inbounds = true function tpic_3d_p2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]*(mpD.v[p,dim]+mpD.∇vᵢⱼ[dim,1,p]*mpD.δnp[nn,1,p]+mpD.∇vᵢⱼ[dim,2,p]*mpD.δnp[nn,2,p]+mpD.∇vᵢⱼ[dim,3,p]*mpD.δnp[nn,3,p])
                if dim == 1
                    @atom meD.mn[no      ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[5,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[4,p])
                elseif dim == 3
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[3,p])
                end
            end
        end
    end
end
function p2n(mpD,meD,g,Δt,instr)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.oobf.= 0.0
    # mapping to mesh
    instr[:cairn][:mapsto][:map].p2n!(ndrange=mpD.nmp,mpD,meD,g);sync(CPU())
    return nothing
end