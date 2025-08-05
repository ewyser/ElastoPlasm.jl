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
@views @kernel inbounds = true function finite_elast(mp::Point{T1,T2},Del) where {T1,T2}
    p = @index(Global)
    if p ≤ mp.nmp 
        # update left cauchy-green tensor
        mp.s.Bᵢⱼ[:,:,p].= mp.s.ΔFᵢⱼ[:,:,p]*mp.s.Bᵢⱼ[:,:,p]*mp.s.ΔFᵢⱼ[:,:,p]'
        # compute logarithmic strain tensor
        λ,n             = eigen(mp.s.Bᵢⱼ[:,:,p],sortby=nothing)
        mp.s.ϵᵢⱼ[:,:,p].= T2(0.5).*(n*diagm(log.(λ))*n')
        # krichhoff stress tensor
        mp.s.τᵢ[:,p]    = Del*mutate(mp.s.ϵᵢⱼ[:,:,p],T2(2.0),:voigt)
    end
end
@views @kernel inbounds = true function infinitesimal_elast(mp::Point{T1,T2},Del) where {T1,T2}
    p = @index(Global)
    if p ≤ mp.nmp 
        # calculate elastic strains & spins
        mp.s.ϵᵢⱼ[:,:,p] .= T2(0.5).*(mp.s.ΔFᵢⱼ[:,:,p]+mp.s.ΔFᵢⱼ[:,:,p]').-mp.s.δᵢⱼ[:,:] 
        mp.s.ωᵢⱼ[:,:,p] .= T2(0.5).*(mp.s.ΔFᵢⱼ[:,:,p]-mp.s.ΔFᵢⱼ[:,:,p]')
        # update cauchy stress tensor
        mp.s.σJᵢⱼ[:,:,p].= mutate(mp.s.σᵢ[:,p],T2(1.0),:tensor)
        mp.s.σJᵢⱼ[:,:,p].= mp.s.σJᵢⱼ[:,:,p]*mp.s.ωᵢⱼ[:,:,p]'+mp.s.σJᵢⱼ[:,:,p]'*mp.s.ωᵢⱼ[:,:,p]
        mp.s.σᵢ[:,p]   .+= Del*mutate(mp.s.ϵᵢⱼ[:,:,p],T2(2.0),:voigt).+mutate(mp.s.σJᵢⱼ[:,:,p],T2(1.0),:voigt)
    end  
end

function init_elast(instr::Dict)
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


function elast(mp::Point{T1,T2},cmpr::NamedTuple,instr::Dict) where {T1,T2}
    instr[:cairn][:elastoplast][:elast].elast!(ndrange=mp.nmp,mp,cmpr.Del);sync(CPU())
    return nothing
end



















































@views @kernel inbounds = true function ELAST(mp,Del,instr)
    p = @index(Global)
    # deformation framework dispatcher
    if instr[:fwrk] == "finite"
        if p ≤ mp.nmp 
            # update left cauchy-green tensor
            mp.s.Bᵢⱼ[:,:,p].= mp.s.ΔFᵢⱼ[:,:,p]*mp.s.Bᵢⱼ[:,:,p]*mp.s.ΔFᵢⱼ[:,:,p]'
            # compute logarithmic strain tensor
            λ,n             = eigen(mp.s.Bᵢⱼ[:,:,p],sortby=nothing)
            mp.s.ϵᵢⱼ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mp.s.τᵢ[:,p]    = Del*mutate(mp.s.ϵᵢⱼ[:,:,p],2.0,:voigt)
        end
    elseif instr[:fwrk] == "infinitesimal"
        if p ≤ mp.nmp 
            # calculate elastic strains & spins
            for i ∈ 1:mp.ndim , j ∈ 1:mp.ndim
                mp.s.ϵᵢⱼ[i,j,p] = 0.5*(mp.s.ΔFᵢⱼ[i,j,p]+mp.s.ΔFᵢⱼ[j,i,p])-mp.s.δᵢⱼ[i,j] 
                mp.s.ωᵢⱼ[i,j,p] = 0.5*(mp.s.ΔFᵢⱼ[i,j,p]-mp.s.ΔFᵢⱼ[j,i,p])
            end
            Δϵ = mutate(mp.s.ϵᵢⱼ[:,:,p],2.0,:voigt)
            Δω = mp.s.ωᵢⱼ[1,2,p]
            σ0 = [Δω*(2.0*mp.s.σᵢ[3,p]),Δω*(-2.0*mp.s.σᵢ[3,p]),Δω*(mp.s.σᵢ[2,p]-mp.s.σᵢ[1,p])]
            for k ∈ 1:length(σ0)
                mp.σᵢ[k,p]+= (Del[k,1]*Δϵ[1]+Del[k,2]*Δϵ[2]+Del[k,3]*Δϵ[3])+σ0[k]
            end
        end   
    end
end