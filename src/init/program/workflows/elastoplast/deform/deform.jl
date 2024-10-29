@views @kernel inbounds = true function measure(mpD,meD,Δt)
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
function init_deformation(instr)
    return measure(CPU())
end
function deform(mpD,meD,Δt,instr)
    instr[:cairn][:elastoplast].deform!(ndrange=mpD.nmp,mpD,meD,Δt);sync(CPU())
    return nothing
end

































