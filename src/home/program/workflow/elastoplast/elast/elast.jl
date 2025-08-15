"""
    mutate(ϵ, Χ, type)

Convert between tensor and Voigt notation for strain or stress, or mutate tensor forms.

# Arguments
- `ϵ`: Strain or stress tensor or vector.
- `Χ`: Scaling factor (e.g., 1.0 or 2.0).
- `type`: Either `:tensor` or `:voigt` for conversion type.

# Returns
- `ϵmut`: Mutated tensor or vector in the requested form.
"""
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
"""
    finite_elast(mpts::Point{T1,T2}, Del) where {T1,T2}

Kernel for finite deformation elasticity update at material points.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `Del`: Elasticity matrix.

# Returns
- Updates stress and strain fields in-place.
"""
@views @kernel inbounds = true function finite_elast(mpts::Point{T1,T2},Del) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # update left cauchy-green tensor
        mpts.s.bᵢⱼ[:,:,p].= mpts.s.ΔFᵢⱼ[:,:,p]*mpts.s.bᵢⱼ[:,:,p]*mpts.s.ΔFᵢⱼ[:,:,p]'
        # compute logarithmic strain tensor
        λ,n             = eigen(mpts.s.bᵢⱼ[:,:,p],sortby=nothing)
        mpts.s.ϵᵢⱼ[:,:,p].= T2(0.5).*(n*diagm(log.(λ))*n')
        # krichhoff stress tensor
        mpts.s.τᵢ[:,p]    = Del*mutate(mpts.s.ϵᵢⱼ[:,:,p],T2(2.0),:voigt)
    end
end
"""
    infinitesimal_elast(mpts::Point{T1,T2}, Del) where {T1,T2}

Kernel for infinitesimal (small strain) elasticity update at material points.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `Del`: Elasticity matrix.

# Returns
- Updates stress and strain fields in-place.
"""
@views @kernel inbounds = true function infinitesimal_elast(mpts::Point{T1,T2},Del) where {T1,T2}
    p = @index(Global)
    if p ≤ mpts.nmp 
        # calculate elastic strains & spins
        mpts.s.ϵᵢⱼ[:,:,p] .= T2(0.5).*(mpts.s.ΔFᵢⱼ[:,:,p]+mpts.s.ΔFᵢⱼ[:,:,p]').-mpts.δᵢⱼ[:,:] 
        mpts.s.ωᵢⱼ[:,:,p] .= T2(0.5).*(mpts.s.ΔFᵢⱼ[:,:,p]-mpts.s.ΔFᵢⱼ[:,:,p]')
        # update cauchy stress tensor
        mpts.s.σJᵢⱼ[:,:,p].= mutate(mpts.s.σᵢ[:,p],T2(1.0),:tensor)
        mpts.s.σJᵢⱼ[:,:,p].= mpts.s.σJᵢⱼ[:,:,p]*mpts.s.ωᵢⱼ[:,:,p]'+mpts.s.σJᵢⱼ[:,:,p]'*mpts.s.ωᵢⱼ[:,:,p]
        mpts.s.σᵢ[:,p]   .+= Del*mutate(mpts.s.ϵᵢⱼ[:,:,p],T2(2.0),:voigt).+mutate(mpts.s.σJᵢⱼ[:,:,p],T2(1.0),:voigt)
    end  
end

"""
    init_elast(instr::NamedTuple)

Initialize the elasticity kernel based on the deformation framework and performance mode.

# Arguments
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- Named tuple with the selected elasticity kernel.
"""
function init_elast(instr::NamedTuple)
    if instr[:perf][:status]
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


"""
    elast(mpts::Point{T1,T2}, cmpr::NamedTuple, instr::NamedTuple) where {T1,T2}

Dispatch and execute the elasticity kernel for the current simulation step.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function elast(mpts::Point{T1,T2},cmpr::NamedTuple,instr::NamedTuple) where {T1,T2}
    instr[:cairn][:elastoplast][:elast].elast!(ndrange=mpts.nmp,mpts,cmpr.Del);sync(CPU())
    return nothing
end



















































"""
    ELAST(mpts, Del, instr)

General elasticity kernel dispatcher for finite or infinitesimal deformation frameworks.

# Arguments
- `mpts`: Material point data structure.
- `Del`: Elasticity matrix.
- `instr`: Instruction/configuration dictionary.

# Returns
- Updates stress and strain fields in-place.
"""
@views @kernel inbounds = true function ELAST(mpts,Del,instr)
    p = @index(Global)
    # deformation framework dispatcher
    if instr[:fwrk] == "finite"
        if p ≤ mpts.nmp 
            # update left cauchy-green tensor
            mpts.s.bᵢⱼ[:,:,p].= mpts.s.ΔFᵢⱼ[:,:,p]*mpts.s.bᵢⱼ[:,:,p]*mpts.s.ΔFᵢⱼ[:,:,p]'
            # compute logarithmic strain tensor
            λ,n             = eigen(mpts.s.bᵢⱼ[:,:,p],sortby=nothing)
            mpts.s.ϵᵢⱼ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mpts.s.τᵢ[:,p]    = Del*mutate(mpts.s.ϵᵢⱼ[:,:,p],2.0,:voigt)
        end
    elseif instr[:fwrk] == "infinitesimal"
        if p ≤ mpts.nmp 
            # calculate elastic strains & spins
            for i ∈ 1:mpts.ndim , j ∈ 1:mpts.ndim
                mpts.s.ϵᵢⱼ[i,j,p] = 0.5*(mpts.s.ΔFᵢⱼ[i,j,p]+mpts.s.ΔFᵢⱼ[j,i,p])-mpts.s.δᵢⱼ[i,j] 
                mpts.s.ωᵢⱼ[i,j,p] = 0.5*(mpts.s.ΔFᵢⱼ[i,j,p]-mpts.s.ΔFᵢⱼ[j,i,p])
            end
            Δϵ = mutate(mpts.s.ϵᵢⱼ[:,:,p],2.0,:voigt)
            Δω = mpts.s.ωᵢⱼ[1,2,p]
            σ0 = [Δω*(2.0*mpts.s.σᵢ[3,p]),Δω*(-2.0*mpts.s.σᵢ[3,p]),Δω*(mpts.s.σᵢ[2,p]-mpts.s.σᵢ[1,p])]
            for k ∈ 1:length(σ0)
                mpts.σᵢ[k,p]+= (Del[k,1]*Δϵ[1]+Del[k,2]*Δϵ[2]+Del[k,3]*Δϵ[3])+σ0[k]
            end
        end   
    end
end