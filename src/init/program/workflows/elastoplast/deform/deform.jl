@views @kernel inbounds = true function MEASURE(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp 
        # compute velocity & displacement gradients
        mpD.∇vᵢⱼ[:,:,p].= 0.0
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            for i ∈ 1:meD.nD , j ∈ 1:meD.nD
                mpD.∇vᵢⱼ[i,j,p]+= mpD.ϕ∂ϕ[nn,p,j+1]*meD.vn[no,i]
            end
        end
        # compute incremental deformation and update
        mpD.ΔFᵢⱼ[:,:,p].= mpD.δᵢⱼ+(Δt.*mpD.∇vᵢⱼ[:,:,p])
        mpD.Fᵢⱼ[:,:,p] .= mpD.ΔFᵢⱼ[:,:,p]*mpD.Fᵢⱼ[:,:,p]
        # update material point's volume
        mpD.ΔJ[p]       = det(mpD.ΔFᵢⱼ[:,:,p])
        mpD.J[p]        = det(mpD.Fᵢⱼ[:,:,p])
        mpD.Ω[p]        = mpD.J[p]*mpD.Ω₀[p]
    end
end
@views @kernel inbounds = true function measure(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp 
        # compute velocity & displacement gradients
        nn              = findall(x->x!=0,mpD.p2n[:,p])
        mpD.∇vᵢⱼ[:,:,p].= (permutedims(mpD.ϕ∂ϕ[nn,p,2:end],(2,1))*meD.vn[meD.e2n[nn,mpD.p2e[p]],:])'
        mpD.∇uᵢⱼ[:,:,p].= Δt.*mpD.∇vᵢⱼ[:,:,p]
        # compute incremental deformation gradient
        mpD.ΔFᵢⱼ[:,:,p].= mpD.δᵢⱼ.+mpD.∇uᵢⱼ[:,:,p]
        mpD.ΔJ[p]       = det(mpD.ΔFᵢⱼ[:,:,p])
        # update deformation gradient
        mpD.Fᵢⱼ[:,:,p] .= mpD.ΔFᵢⱼ[:,:,p]*mpD.Fᵢⱼ[:,:,p]
        # update material point's volume
        mpD.J[p]        = det(mpD.Fᵢⱼ[:,:,p])
        mpD.Ω[p]        = mpD.J[p]*mpD.Ω₀[p]
    end
end
function init_deformation(instr)
    if instr[:perf]
        return MEASURE(CPU())
    else
        return measure(CPU())
    end
    return nothing
end
function deformation!(mpD,meD,Δt,instr)
    instr[:cairn][:elastoplast].deform!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
    return nothing
end

































