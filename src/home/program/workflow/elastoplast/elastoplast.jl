function elasto(mp::Point{T1,T2},mesh,cmpr::NamedTuple,dt::T2,instr::Dict) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    update(mp,mesh,dt,instr)
    # update {kirchoff|cauchy} stresses
    elast(mp,cmpr,instr)
    return nothing
end
function elastoplast(mp::Point{T1,T2},mesh,cmpr::NamedTuple,dt::T2,instr::Dict) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    update(mp,mesh,dt,instr)
    # update {kirchoff|cauchy} stresses
    elast(mp,cmpr,instr)
    # plastic corrector
    ηmax = plast(mp,mesh,cmpr,instr)
    return ηmax::Int64
end