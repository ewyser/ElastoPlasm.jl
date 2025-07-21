@views function mutate(ϵ,Χ,type)
    if type == :tensor # α = 1/2 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (3,)
            ϵmut = [  ϵ[1] Χ*ϵ[3];
                    Χ*ϵ[3]   ϵ[2]]
        elseif size(ϵ) == (6,)
            ϵmut = [  ϵ[1] Χ*ϵ[6] Χ*ϵ[5];
                    Χ*ϵ[6]   ϵ[2] Χ*ϵ[4];
                    Χ*ϵ[5] Χ*ϵ[4]   ϵ[3]]
        end
    elseif type == :voigt # α = 2.0 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (2,2)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],Χ*ϵ[1,2]) #xx,yy,zz,xy
        elseif size(ϵ) == (3,3)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],ϵ[3,3],Χ*ϵ[2,3],Χ*ϵ[1,3],Χ*ϵ[1,2]) #xx,yy,zz,yz,xz,xy
        end
    end
    return ϵmut
end
@views @kernel inbounds = true function ELAST(mp,Del,instr)
    p = @index(Global)
    # deformation framework dispatcher
    if instr[:fwrk] == "finite"
        if p ≤ mp.nmp 
            # update left cauchy-green tensor
            mp.Bᵢⱼ[:,:,p].= mp.ΔFᵢⱼ[:,:,p]*mp.Bᵢⱼ[:,:,p]*mp.ΔFᵢⱼ[:,:,p]'
            # compute logarithmic strain tensor
            λ,n            = eigen(mp.Bᵢⱼ[:,:,p],sortby=nothing)
            mp.ϵᵢⱼ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mp.τᵢ[:,p]    = Del*mutate(mp.ϵᵢⱼ[:,:,p],2.0,:voigt)
        end
    elseif instr[:fwrk] == "infinitesimal"
        if p ≤ mp.nmp 
            # calculate elastic strains & spins
            for i ∈ 1:mp.ndim , j ∈ 1:mp.ndim
                mp.ϵᵢⱼ[i,j,p] = 0.5*(mp.ΔFᵢⱼ[i,j,p]+mp.ΔFᵢⱼ[j,i,p])-mp.δᵢⱼ[i,j] 
                mp.ωᵢⱼ[i,j,p] = 0.5*(mp.ΔFᵢⱼ[i,j,p]-mp.ΔFᵢⱼ[j,i,p])
            end
            Δϵ = mutate(mp.ϵᵢⱼ[:,:,p],2.0,:voigt)
            Δω = mp.ωᵢⱼ[1,2,p]
            σ0 = [Δω*(2.0*mp.σᵢ[3,p]),Δω*(-2.0*mp.σᵢ[3,p]),Δω*(mp.σᵢ[2,p]-mp.σᵢ[1,p])]
            for k ∈ 1:length(σ0)
                mp.σᵢ[k,p]+= (Del[k,1]*Δϵ[1]+Del[k,2]*Δϵ[2]+Del[k,3]*Δϵ[3])+σ0[k]
            end
        end   
    end
end
@views @kernel inbounds = true function finite_elast(mp,Del,instr)
    p = @index(Global)
    if p ≤ mp.nmp 
        # update left cauchy-green tensor
        mp.Bᵢⱼ[:,:,p].= mp.ΔFᵢⱼ[:,:,p]*mp.Bᵢⱼ[:,:,p]*mp.ΔFᵢⱼ[:,:,p]'
        # compute logarithmic strain tensor
        λ,n            = eigen(mp.Bᵢⱼ[:,:,p],sortby=nothing)
        mp.ϵᵢⱼ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
        # krichhoff stress tensor
        mp.τᵢ[:,p]    = Del*mutate(mp.ϵᵢⱼ[:,:,p],2.0,:voigt)
    end
end
@views @kernel inbounds = true function infinitesimal_elast(mp,Del,instr)
    p = @index(Global)
    if p ≤ mp.nmp 
        # calculate elastic strains & spins
        mp.ϵᵢⱼ[:,:,p] .= 0.5.*(mp.ΔFᵢⱼ[:,:,p]+mp.ΔFᵢⱼ[:,:,p]').-mp.δᵢⱼ[:,:] 
        mp.ωᵢⱼ[:,:,p] .= 0.5.*(mp.ΔFᵢⱼ[:,:,p]-mp.ΔFᵢⱼ[:,:,p]')
        # update cauchy stress tensor
        mp.σJᵢⱼ[:,:,p].= mutate(mp.σᵢ[:,p],1.0,:tensor)
        mp.σJᵢⱼ[:,:,p].= mp.σJᵢⱼ[:,:,p]*mp.ωᵢⱼ[:,:,p]'+mp.σJᵢⱼ[:,:,p]'*mp.ωᵢⱼ[:,:,p]
        mp.σᵢ[:,p]   .+= Del*mutate(mp.ϵᵢⱼ[:,:,p],2.0,:voigt).+mutate(mp.σJᵢⱼ[:,:,p],1.0,:voigt)
    end  
end

function init_elast(instr)
    if instr[:perf]
        kernel1 = ELAST(CPU())
    else
        # deformation framework dispatcher
        if instr[:fwrk][:deform] == "finite"
            kernel1 = finite_elast(CPU())
        elseif instr[:fwrk][:deform] == "infinitesimal"
            kernel1 = infinitesimal_elast(CPU())
        end
    end
    return (;elast! = kernel1,)
end


function elast(mp,cmParam,instr,type)
    instr[:cairn][:elastoplast][:elast].elast!(ndrange=mp.nmp,mp,cmParam.Del,instr);sync(CPU())
    return nothing
end