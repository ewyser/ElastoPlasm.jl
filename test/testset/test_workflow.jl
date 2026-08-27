"""
    generate_geometry_cases() -> Vector{NamedTuple}

Generate all combinations of shape functions.
"""
function generate_geometry_cases()
    return [
        (dim = 2, L = [64.1584, 64.1584 / 4.0]               , nel = [40, 10]    ),
        #(dim = 3, L = [64.1584, 64.1584 / 4.0, 64.1584 / 4.0], nel = [40, 10, 10]),
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
    generate_strain_cases() -> Vector{NamedTuple}

Generate all `strain` config cases (deformation framework).
"""
function generate_strain_cases()
    return [(deform = d,) for d in ["finite", "infinitesimal"]]
end

"""
    generate_transfer_cases() -> Vector{NamedTuple}

Generate all `basis.trsfr`/`basis.C_pf` cases (P2G/G2P transfer scheme).
"""
function generate_transfer_cases()
    return [(trsfr = t, C_pf = 1.0) for t in ["std", "tpic", "apic"]]
end

"""
    generate_stab_cases() -> Vector{NamedTuple}

Generate all `stab` config cases (locking correction, MUSL reprojection).
"""
function generate_stab_cases()
    return [(locking = l, damping = 0.1, musl = m) for l in [true, false] for m in [true, false]]
end

"""
    generate_nonloc_cases() -> Vector{NamedTuple}

Generate all `nonloc` config cases (nonlocal plastic-strain regularization).
"""
function generate_nonloc_cases()
    return [(status = s, ls = 0.5) for s in [true, false]]
end

@testset "+ $(basename(@__FILE__))" verbose = true begin

    cfg_geom     = generate_geometry_cases()
    cfg_basis    = generate_basis_cases()
    cfg_strain   = generate_strain_cases()
    cfg_transfer = generate_transfer_cases()
    cfg_stab     = generate_stab_cases()
    cfg_nonloc   = generate_nonloc_cases()
    cfg_grf      = (;
        status = false,
        covariance = "gaussian",
        param = (;Iₓ = [2.5, 2.5, 2.5], Nₕ = 5000, kₘ = 100,),
    )

    sim  = 1
    nsim = length(cfg_geom) * length(cfg_basis) * length(cfg_strain) * length(cfg_transfer) * length(cfg_stab) * length(cfg_nonloc)
    prog = ElastoPlasm.make_progress(nsim; desc="Executing $(basename(@__FILE__))")
    for (k, geom) in enumerate(cfg_geom)
        @testset "$(geom.dim)d geometry" verbose = true begin
            for (l, basis) in enumerate(cfg_basis)
                @testset "$(basis.which) basis" verbose = true begin
                    for strain in cfg_strain, transfer in cfg_transfer, stab in cfg_stab, nonloc in cfg_nonloc
                        @testset "$(strain.deform), $(transfer.trsfr), locking=$(stab.locking), musl=$(stab.musl), nonloc=$(nonloc.status)" verbose = true begin
                            name = "slump_$(geom.dim)d_$(basis.which)_$(strain.deform)_$(transfer.trsfr)_lock$(stab.locking)_musl$(stab.musl)_nl$(nonloc.status)"
                            status = false
                            try
                                @suppress begin
                                    basiscfg = merge(basis, transfer)
                                    jld2     = slump_problem(geom.L, geom.nel; fid="test/$(name)", grf=cfg_grf, basis=basiscfg, strain=strain, stab=stab, nonloc=nonloc)
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