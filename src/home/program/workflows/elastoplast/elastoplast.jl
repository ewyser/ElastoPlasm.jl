function elastoplast(mp,mesh,cmpr,dt,instr)
    # update {logarithmic|infinitesimal} strains
    update(mp,mesh,dt,instr)
    # update {kirchoff|cauchy} stresses
    elast(mp,cmpr,instr)
    # plastic corrector
    ηmax = plast(mp,mesh,cmpr,instr)
    return ηmax::Int64
end