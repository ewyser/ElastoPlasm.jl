function init_update(instr::NamedTuple)
    if instr[:fwrk][:deform] == "finite"
        kernel1 = finite_deform(CPU())
    elseif instr[:fwrk][:deform] == "infinitesimal"
        kernel1 = infinitesimal_deform(CPU())
    end
    if instr[:basis][:which] == "gimpm"
        if instr[:fwrk][:deform] == "finite"
            if instr[:basis][:how] == "undeformed"
                kernel2 = undeformed(CPU())
            elseif instr[:basis][:how] == "detFij"
                kernel2 = detFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "Fii"
                kernel2 = Fᵢᵢ(CPU())
            elseif instr[:basis][:how] == "Uii"
                kernel2 = Uᵢᵢ(CPU())
            end
        elseif instr[:fwrk][:deform] == "infinitesimal"
            if instr[:basis][:how] == "undeformed"
                kernel2 = undeformed(CPU())
            elseif instr[:basis][:how] == "detΔFij" 
                kernel2 = detΔFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "ΔFii"
                kernel2 = ΔFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "ΔUii"
                kernel2 = ΔUᵢᵢ(CPU())
            end
        end
    else
        kernel2 = nothing
    end
    if instr[:fwrk][:locking]
        kernel3a = ΔJn(CPU())
        kernel3b = ΔJs(CPU())
        kernel3c = ΔJp(CPU())
    else
        kernel3a = nothing
        kernel3b = nothing
        kernel3c = nothing
    end
    return (;deform! = kernel1,domain! = kernel2,ΔJn! = kernel3a,ΔJs! = kernel3b,ΔJp! = kernel3c,)
end
function update(mpts::Point{T1,T2},mesh::Mesh{T1,T2},dt::T2,instr::NamedTuple) where {T1,T2}
    # get incremental deformation tensor
    instr[:cairn][:elastoplast][:update].deform!(ndrange=mpts.nmp,mpts,mesh,dt);sync(CPU())
    # update material point's domain
    if instr[:basis][:which] == "gimpm"
        instr[:cairn][:elastoplast][:update].domain!(ndrange=mpts.nmp,mpts);sync(CPU())
    end
    # volumetric locking correction
    if instr[:fwrk][:locking]
        # init mesh quantities to zero
        mesh.ΔJ.= T2(0.0)
        # calculate dimensional cst.
        dim     = T2(1.0)/mesh.dim
        # mapping to mesh 
        instr[:cairn][:elastoplast][:update].ΔJn!(ndrange=mpts.nmp,mpts,mesh);sync(CPU())
        # compute nodal determinant of incremental deformation 
        instr[:cairn][:elastoplast][:update].ΔJs!(ndrange=mesh.nno[end],mesh);sync(CPU())
        # compute determinant Jbar 
        instr[:cairn][:elastoplast][:update].ΔJp!(ndrange=mpts.nmp,mpts,mesh,dim);sync(CPU())
    end  
    return nothing
end
