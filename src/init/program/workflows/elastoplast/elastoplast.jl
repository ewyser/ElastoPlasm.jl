function elastoplast(mp,mesh,cmParam,Δt,instr)
    # update {logarithmic|infinitesimal} strains
    update(mp,mesh,Δt,instr)
    # update {kirchoff|cauchy} stresses
    elast(mp,cmParam,instr,:update)
    # plastic corrector
    ηmax = plast(mp,mesh,cmParam,instr)
    return ηmax::Int64
end