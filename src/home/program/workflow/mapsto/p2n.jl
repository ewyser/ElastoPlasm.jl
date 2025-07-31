@views @kernel inbounds = true function transform(mp)
    p = @index(Global)
    # deformation framework dispatcher
    if p ≤ mp.nmp 
        mp.σᵢ[:,p].= mp.τᵢ[:,p]./mp.J[p]
    end   
end
@kernel inbounds = true function flip_1d_p2n(mp,mesh,g)
    p = @index(Global)
    if p≤mp.nmp 
        # accumulation
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            @atom mesh.p[no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*mp.v[p])
            # lumped mass matrix
            @atom mesh.mᵢ[no]+= mp.ϕ∂ϕ[nn,p,1]*mp.m[p]
            # consistent mass matrix
            # mesh.Mᵢⱼ[mp.p2n[:,p],mp.p2n[:,p]].+= (mp.ϕ∂ϕ[:,p,1].*mp.ϕ∂ϕ[:,p,1]').*mp.m[p]   
            @atom mesh.oobf[no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*g) 
            @atom mesh.oobf[no]-= mp.Ω[p]*mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[1,p]
        end
    end
end
@kernel inbounds = true function flip_2d_p2n(mp,mesh,g)
    p = @index(Global)
    if p≤mp.nmp
        for dim ∈ 1:mesh.dim 
            # accumulation
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if no < 1 continue end
                @atom mesh.p[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*mp.v[dim,p])
                if dim == 1
                    # lumped mass matrix
                    @atom mesh.mᵢ[no]+= mp.ϕ∂ϕ[nn,p,1]*mp.m[p]
                    # consistent mass matrix
                    # mesh.Mᵢⱼ[mp.p2n[:,p],mp.p2n[:,p]].+= (mp.ϕ∂ϕ[:,p,1].*mp.ϕ∂ϕ[:,p,1]').*mp.m[p]    
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[1,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[3,p])
                elseif dim == 2
                    @atom mesh.oobf[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*g[dim]      )
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[3,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[2,p])
                end
            end
        end
    end
end
@kernel inbounds = true function flip_3d_p2n(mp,mesh,g)
    p = @index(Global)
    if p≤mp.nmp
        for dim ∈ 1:mesh.dim 
            # accumulation
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if no < 1 continue end
                @atom mesh.p[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*mp.v[dim,p])
                if dim == 1
                    @atom mesh.mᵢ[no      ]+= mp.ϕ∂ϕ[nn,p,1]*mp.m[p] 
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[1,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[6,p]+mp.ϕ∂ϕ[nn,p,4]*mp.σᵢ[5,p])
                elseif dim == 2
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[6,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[2,p]+mp.ϕ∂ϕ[nn,p,4]*mp.σᵢ[4,p])
                elseif dim == 3
                    @atom mesh.oobf[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*g[dim]      )
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[5,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[4,p]+mp.ϕ∂ϕ[nn,p,4]*mp.σᵢ[3,p])
                end
            end
        end
    end
end
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TPIC transfer scheme, see Nakamura etal, 2023, https://doi.org/10.1016/j.cma.2022.115720
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
@kernel inbounds = true function tpic_2d_p2n(mp,mesh,g)
    p = @index(Global)
    if p≤mp.nmp
        for dim ∈ 1:mesh.dim 
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if no < 1 continue end
                @atom mesh.p[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*mp.m[p]*(mp.v[dim,p]+mp.∇vᵢⱼ[dim,1,p]*mp.δnp[nn,1,p]+mp.∇vᵢⱼ[dim,2,p]*mp.δnp[nn,2,p])
                if dim == 1
                    @atom mesh.mᵢ[no]      += mp.ϕ∂ϕ[nn,p,1]*mp.m[p]
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[1,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[3,p])
                elseif dim == 2
                    @atom mesh.oobf[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*g[dim]      )
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[3,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[2,p])
                end
            end
        end
    end

end
@kernel inbounds = true function tpic_3d_p2n(mp,mesh,g)
    p = @index(Global)
    if p≤mp.nmp
        for dim ∈ 1:mesh.dim 
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if no < 1 continue end
                @atom mesh.p[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*mp.m[p]*(mp.v[dim,p]+mp.∇vᵢⱼ[dim,1,p]*mp.δnp[nn,1,p]+mp.∇vᵢⱼ[dim,2,p]*mp.δnp[nn,2,p]+mp.∇vᵢⱼ[dim,3,p]*mp.δnp[nn,3,p])
                if dim == 1
                    @atom mesh.mᵢ[no      ]+= mp.ϕ∂ϕ[nn,p,1]*mp.m[p]
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[1,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[6,p]+mp.ϕ∂ϕ[nn,p,4]*mp.σᵢ[5,p])
                elseif dim == 2
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[6,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[2,p]+mp.ϕ∂ϕ[nn,p,4]*mp.σᵢ[4,p])
                elseif dim == 3
                    @atom mesh.oobf[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*g[dim]      )
                    @atom mesh.oobf[dim,no]-= mp.Ω[p]*(mp.ϕ∂ϕ[nn,p,2]*mp.σᵢ[5,p]+mp.ϕ∂ϕ[nn,p,3]*mp.σᵢ[4,p]+mp.ϕ∂ϕ[nn,p,4]*mp.σᵢ[3,p])
                end
            end
        end
    end
end
function p2n(mp,mesh,g,dt,instr)
    # get cauchy stress 
    if instr[:fwrk][:deform] == "finite"
        instr[:cairn][:mapsto][:map].σᵢ!(ndrange=mp.nmp,mp);sync(CPU())
    end
    # initialize nodal quantities
    mesh.mᵢ  .= 0.0
    mesh.p  .= 0.0
    mesh.oobf.= 0.0
    # mapping to mesh
    instr[:cairn][:mapsto][:map].p2n!(ndrange=mp.nmp,mp,mesh,g);sync(CPU())
    return nothing
end