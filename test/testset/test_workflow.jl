"""
    generate_geometry_cases() -> Vector{NamedTuple}

Generate all combinations of shape functions.
"""
function generate_geometry_cases()
    return [
        (dim = 2, L = [64.1584, 64.1584 / 4.0]               , nel = [40, 10]    ),
        (dim = 3, L = [64.1584, 64.1584 / 4.0, 64.1584 / 4.0], nel = [40, 10, 10]),
    ]
end


"""
    generate_basis_cases() -> Vector{NamedTuple}

Generate all combinations of shape functions.
"""
function generate_basis_cases()
    return [
        (which = "smpm"  , how = nothing     ),
        (which = "gimpm" , how = "undeformed"),
        (which = "bsmpm" , how = nothing     ),
        (which = "mlsmpm", how = nothing     ),
    ]
end

"""
    generate_fwrk_cases() -> Vector{NamedTuple}

Generate all combinations of deformation frameworks, transfer schemes, and locking options.
"""
function generate_fwrk_cases()
    return [(deform = d, trsfr = t, locking = l, C_pf = 1.0, musl = m, damping = 0.1) 
            for d in ["finite", "infinitesimal"] 
            for t in ["std", "tpic", "apic"] 
            for l in [true, false]
            for m in [true, false]
        ]
end

@testset "+ $(basename(@__FILE__))" verbose = true begin

    cfg_geom  = generate_geometry_cases()
    cfg_fwrk  = generate_fwrk_cases()
    cfg_basis = generate_basis_cases()
    cfg_grf   = (;
        status = false,
        covariance = "gaussian",
        param = (;Iₓ = [2.5, 2.5, 2.5], Nₕ = 5000, kₘ = 100,),
    )
    cfg_nonloc = (;
        status = false,
        ls = 0.5,
    )

    sim  = 1
    nsim = length(cfg_geom) * length(cfg_fwrk) * length(cfg_basis)
    prog = ElastoPlasm.make_progress(nsim; desc="Executing $(basename(@__FILE__))")
    for (k, geom) in enumerate(cfg_geom)
        @testset "$(geom.dim)d geometry" verbose = true begin 
            for (l, basis) in enumerate(cfg_basis)
                @testset "$(basis.which) basis" verbose = true begin
                    for (m, fwrk) in enumerate(cfg_fwrk)
                        @testset "$(fwrk.deform), $(fwrk.trsfr), locking=$(fwrk.locking), musl=$(fwrk.musl)" verbose = true begin
                            name = "slump_$(geom.dim)d_$(basis.which)_$(fwrk.deform)_$(fwrk.trsfr)_lock$(fwrk.locking)_musl$(fwrk.musl)"
                            status = false
                            try
                                @suppress begin
                                    strain   = (; deform = fwrk.deform)
                                    basiscfg = merge(basis, (; trsfr = fwrk.trsfr, C_pf = fwrk.C_pf))
                                    stab     = (; locking = fwrk.locking, damping = fwrk.damping, musl = fwrk.musl)
                                    jld2     = slump_problem(geom.L, geom.nel; fid="test/$(name)", grf=cfg_grf, basis=basiscfg, strain=strain, stab=stab, nonloc=cfg_nonloc)
                                    status   = elastoplasm!(jld2; workflows=[elastodynamic!, elastoplastic!]).success
                                end
                            catch e
                                @warn "Test case failed" name exception=(e, catch_backtrace())
                            end
                            @test status
                            sim += 1
                            next!(prog; desc="Executing $(basename(@__FILE__)) $(sim)/$nsim...")
                        end
                    end  
                end   
            end
        end
    end
    #metanalysis(joinpath(ElastoPlasm.self.sys.dump,"test"))
end