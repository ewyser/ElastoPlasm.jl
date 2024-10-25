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
@views @kernel inbounds = true function ELAST!(mpD,Del,instr)
    p = @index(Global)
    # deformation framework dispatcher
    if instr[:fwrk] == "finite"
        if p ≤ mpD.nmp 
            # update left cauchy-green tensor
            mpD.Bᵢⱼ[:,:,p].= mpD.ΔFᵢⱼ[:,:,p]*mpD.Bᵢⱼ[:,:,p]*mpD.ΔFᵢⱼ[:,:,p]'
            # compute logarithmic strain tensor
            λ,n            = eigen(mpD.Bᵢⱼ[:,:,p],sortby=nothing)
            mpD.ϵᵢⱼ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mpD.τᵢ[:,p]    = Del*mutate(mpD.ϵᵢⱼ[:,:,p],2.0,:voigt)
        end
    elseif instr[:fwrk] == "infinitesimal"
        if p ≤ mpD.nmp 
            # calculate elastic strains & spins
            for i ∈ 1:mpD.nD , j ∈ 1:mpD.nD
                mpD.ϵᵢⱼ[i,j,p] = 0.5*(mpD.ΔFᵢⱼ[i,j,p]+mpD.ΔFᵢⱼ[j,i,p])-mpD.δᵢⱼ[i,j] 
                mpD.ωᵢⱼ[i,j,p] = 0.5*(mpD.ΔFᵢⱼ[i,j,p]-mpD.ΔFᵢⱼ[j,i,p])
            end
            Δϵ = mutate(mpD.ϵᵢⱼ[:,:,p],2.0,:voigt)
            Δω = mpD.ωᵢⱼ[1,2,p]
            σ0 = [Δω*(2.0*mpD.σᵢ[3,p]),Δω*(-2.0*mpD.σᵢ[3,p]),Δω*(mpD.σᵢ[2,p]-mpD.σᵢ[1,p])]
            for k ∈ 1:length(σ0)
                mpD.σᵢ[k,p]+= (Del[k,1]*Δϵ[1]+Del[k,2]*Δϵ[2]+Del[k,3]*Δϵ[3])+σ0[k]
            end
        end   
    end
end
@views @kernel inbounds = true function elast!(mpD,Del,instr)
    p = @index(Global)
    # deformation framework dispatcher
    if instr[:fwrk] == "finite"
        if p ≤ mpD.nmp 
            # update left cauchy-green tensor
            mpD.Bᵢⱼ[:,:,p].= mpD.ΔFᵢⱼ[:,:,p]*mpD.Bᵢⱼ[:,:,p]*mpD.ΔFᵢⱼ[:,:,p]'
            # compute logarithmic strain tensor
            λ,n            = eigen(mpD.Bᵢⱼ[:,:,p],sortby=nothing)
            mpD.ϵᵢⱼ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mpD.τᵢ[:,p]    = Del*mutate(mpD.ϵᵢⱼ[:,:,p],2.0,:voigt)
        end
    elseif instr[:fwrk] == "infinitesimal"
        if p ≤ mpD.nmp 
            # calculate elastic strains & spins
            mpD.ϵᵢⱼ[:,:,p] .= 0.5.*(mpD.ΔFᵢⱼ[:,:,p]+mpD.ΔFᵢⱼ[:,:,p]').-mpD.δᵢⱼ[i,j] 
            mpD.ωᵢⱼ[:,:,p] .= 0.5.*(mpD.ΔFᵢⱼ[:,:,p]-mpD.ΔFᵢⱼ[:,:,p]')
            # update cauchy stress tensor
            mpD.σJᵢⱼ[:,:,p].= mutate(mpD.σᵢ[:,p],1.0,:tensor)
            mpD.σJᵢⱼ[:,:,p].= mpD.σJᵢⱼ[:,:,p]*mpD.ωᵢⱼ[:,:,p]'+mpD.σJᵢⱼ[:,:,p]'*mpD.ωᵢⱼ[:,:,p]
            mpD.σᵢ[:,p]   .+= Del*mutate(mpD.ϵᵢⱼ[:,:,p],2.0,:voigt).+mutate(mpD.σJᵢⱼ[:,:,p],1.0,:voigt)
        end   
    end
end
@views @kernel inbounds = true function transform!(mpD)
    p = @index(Global)
    # deformation framework dispatcher
    if p ≤ mpD.nmp 
        mpD.σᵢ[:,p] .= mpD.τᵢ[:,p]./mpD.J[p]
    end   
end

function init_stress(instr;what="τorσ!")
    if instr[:perf]
        kernel1 = ELAST!(CPU())
    else
        kernel1 = elast!(CPU())
    end
    kernel2 = transform!(CPU())
    return (;τorσ! = kernel1, cauchy! = kernel2,)
end
function stress(mpD,cmParam,instr,type)
    if type == :update
        instr[:cairn][:elastoplast][:stress].τorσ!(mpD,cmParam.Del,instr; ndrange=mpD.nmp);sync(CPU())
    elseif type == :cauchy
        instr[:cairn][:elastoplast][:stress].cauchy!(mpD; ndrange=mpD.nmp);sync(CPU())
    end
    return nothing
end