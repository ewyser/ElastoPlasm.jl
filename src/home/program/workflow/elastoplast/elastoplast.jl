function elasto(mpts::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::Dict) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    update(mpts,mesh,dt,instr)
    # update {kirchoff|cauchy} stresses
    elast(mpts,cmpr,instr)
    return nothing
end
function elastoplast(mpts::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::Dict) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    update(mpts,mesh,dt,instr)
    # update {kirchoff|cauchy} stresses
    elast(mpts,cmpr,instr)
    # plastic corrector
    plast(mpts,mesh,cmpr,instr)
    return nothing
end