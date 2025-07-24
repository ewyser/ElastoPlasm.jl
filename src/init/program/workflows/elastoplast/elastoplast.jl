function elastoplast(mp,mesh,cmp,dt,instr)
    # update {logarithmic|infinitesimal} strains
    update(mp,mesh,dt,instr)
    # update {kirchoff|cauchy} stresses
    elast(mp,cmp,instr,:update)
    # plastic corrector
    ηmax = plast(mp,mesh,cmp,instr)
    return ηmax::Int64
end