"""
    init_plast(instr)

Initialize plasticity and nonlocal kernels based on the instruction set.

# Arguments
- `instr`: Instruction/configuration dictionary.

# Returns
- Named tuple with nonlocal and return-mapping kernel functions.
"""
function init_plast(instr)
    kernel1 = nonlocal(CPU())
    if instr[:plast][:constitutive] == "MC"
        #ηmax = MCRetMap!(mpts,ϵpII,cmp,instr[:fwrk])
    elseif instr[:plast][:constitutive] == "DP"        
        if instr[:fwrk][:deform] == "finite"
            kernel2 = finite_DP(CPU())
        elseif instr[:fwrk][:deform] == "infinitesimal"
            kernel2 = infinitesimal_DP(CPU())
        end
    elseif instr[:plast][:constitutive] == "J2"
        if instr[:fwrk][:deform] == "finite"
            kernel2 = finite_J2(CPU())
        elseif instr[:fwrk][:deform] == "infinitesimal"
            kernel2 = infinitesimal_J2(CPU())
        end
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mpts,cmp,instr[:fwrk])
    else
        throw(error("InvalidReturnMapping: $(instr[:plast][:constitutive])"))
    end 
    return (;nonloc! = kernel1, retmap! = kernel2,) 
end
"""
    plast(mpts::Point{T1,T2}, mesh::Mesh{T1,T2}, cmpr::NamedTuple, instr::NamedTuple) where {T1,T2}

Apply nonlocal regularization and plastic return-mapping for material points.

# Arguments
- `mpts::Point{T1,T2}`: Material point data structure.
- `mesh::Mesh{T1,T2}`: Mesh data structure.
- `cmpr::NamedTuple`: Constitutive model parameters.
- `instr::NamedTuple`: Instruction/configuration dictionary.

# Returns
- `nothing`. Updates fields in-place.
"""
function plast(mpts::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,instr::NamedTuple) where {T1,T2}
    # nonlocal regularization
    if instr[:nonloc][:status]
        ls             = T2(instr[:nonloc][:ls])
        mpts.e2p        .= T1(0)
        mpts.p2p        .= T1(0)
        mpts.s.ϵpII[2,:].= T2(0.0)
        W,w     = spzeros(T2,mpts.nmp),spzeros(T2,mpts.nmp,mpts.nmp)
        for proc ∈ ["tplgy","p->q","p<-q"]
            instr[:cairn][:elastoplast][:plast].nonloc!(W,w,mpts,mesh,ls,proc; ndrange=mpts.nmp);sync(CPU())
        end
    else
        mpts.s.ϵpII[2,:].= mpts.s.ϵpII[1,:]
    end
    # plastic return-mapping dispatcher
    instr[:cairn][:elastoplast][:plast].retmap!(mpts,cmpr; ndrange=mpts.nmp);sync(CPU())
    return nothing
end