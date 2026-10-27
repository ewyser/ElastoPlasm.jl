@testset "+ $(basename(@__FILE__))" verbose = true begin
    # Every example's initial problem, flattened to plain arrays and compared, bit for bit, with a
    # stored reference (`test/dataset/setup_reference.jld2`, written on the first run): a change under
    # `init/` must leave the set-up particles, constants, boundary conditions and shape functions
    # unchanged. The reference keeps each array's size and a 64-bit FNV-1a digest of its bits, not
    # the arrays. Delete it deliberately when a setup change is intended.
    import Random
    EP     = ElastoPlasm
    noplot = (; plot=(; get_default().plot..., status=false))
    L      = [64.1584, 64.1584/4]
    function flat(jld2)
        pb = load_simulation(jld2).problem; basis, solver = pb.basis, pb.solver
        m, c = pb.mpts, pb.mpts.s.cmp
        return Dict{String,Any}(
            "x" => reduce(hcat, m.x), "Omega" => m.Ω₀, "n0" => m.n₀, "rho0" => m.s.ρ₀, "rho" => m.s.ρ, "m" => m.s.m,
            "Gc" => getfield.(c, :Gc), "Kc" => getfield.(c, :Kc), "Hp" => getfield.(c, :Hp), "c0" => getfield.(c, :c₀),
            "cr" => getfield.(c, :cᵣ), "phi" => hasfield(eltype(c), :ϕ₀) ? getfield.(c, :ϕ₀) : Float64[],
            "v" => reduce(hcat, m.s.v), "T" => isnothing(m.t) ? Float64[] : m.t.T,
            "c" => isnothing(m.t) ? Float64[] : m.t.c, "k" => isnothing(m.t) ? Float64[] : m.t.k,
            "bcs" => Float64.(pb.mesh.s.bcs.status), "N" => basis.N,
            "time" => [pb.time.t; pb.time.te; pb.time.tg; pb.time.tep],
            "solver" => string((solver.solution, solver.basis, solver.stab, solver.bcs, solver.material, solver.dtype)))
    end
    # (problem, is the point order meaningful): the collision's disks may come in any order
    cases = [
        ("slump 2D",   () -> setup_simulation(slump_case(L, [40, 10]; cli()..., noplot...))),
        ("slump 3D",   () -> setup_simulation(slump_case([L; L[2]], [40, 10, 10]; cli()..., noplot...))),
        ("slump GRF",  () -> (Random.seed!(1234); setup_simulation(slump_case(L, [40, 10]; cli()..., noplot..., grf=(; get_default().grf..., status=true))))),
        ("slump fine", () -> setup_simulation(slump_case(L, [80, 20]; cli()..., noplot...))),
        ("collapse",   () -> setup_simulation(collapse_case([40, 20]; noplot...))),
        ("block",      () -> setup_simulation(block_case(; noplot...))),
        ("column",     () -> setup_simulation(column_case([1, 20], 0.0, 1.0e4, 80.0, 50.0; noplot...))),
        ("thermal",    () -> setup_simulation(thermal_case(L, [40, 10]; noplot...))),
        ("collision",  () -> setup_simulation(collision_case([1.0, 1.0], [20, 20]; noplot...))),
    ]
    # FNV-1a over the bits of the values (strings: their bytes)
    fnv(words) = foldl((h, w) -> (h ⊻ w) * 0x00000100000001b3, words; init=0xcbf29ce484222325)
    digest(v::AbstractString) = (0, fnv(UInt64.(codeunits(v))))
    digest(v::AbstractArray)  = (size(v), fnv(reinterpret(UInt64, Float64.(vec(v)))))
    file = joinpath(DATASET, "setup_reference.jld2")
    full = Dict{String,Any}(name => flat(f()) for (name, f) ∈ cases)
    # the collision's disks may come in any order: its points (and per-point fields) sorted by position
    function sorted(d)
        p = sortperm(collect(eachcol(d["x"])))
        return Dict(k => (v isa AbstractVector && length(v) == length(p) ? v[p] :
                          v isa AbstractMatrix && k ∈ ("x", "v") ? v[:, p] : v) for (k, v) ∈ d if k != "N")
    end
    now = Dict{String,Any}(name => Dict(k => digest(v) for (k, v) ∈ (name == "collision" ? sorted(d) : d)) for (name, d) ∈ full)
    if !isfile(file)
        jldsave(file; ref=now)
        @warn "Wrote the setup reference $(file); later runs compare against it"
    end
    ref = load(file, "ref")
    @testset "$name" for (name, _) ∈ cases
        a, b = now[name], ref[name]
        @test sort(collect(keys(a))) == sort(collect(keys(b)))
        for k ∈ sort(collect(keys(b)))
            @test (k, a[k]) == (k, b[k])
        end
    end
end
