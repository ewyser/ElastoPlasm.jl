# The explicit workflows on a device, with JLArrays (GPUArrays' reference array type, on the CPU: a KernelAbstractions
# backend whose launches adapt their arguments, and indexing such an array element by element from the host is an
# error), registered as a UnifiedBackend device. It runs the same kernels in the same order as the host, but compiled
# for other argument types, so StaticArrays' `muladd`s fuse into FMAs differently: `F = ΔF·F` in `deform` first differs
# by an ulp, and a short plastic slump by at most 9e-12 relative (Float64), 8e-7 (Float32); integer and Bool arrays
# (connectivity, boundaries) are equal. JLArrays is a test-only dependency (`Pkg.test()` installs it).
const HAS_JLARRAYS = !isnothing(Base.find_package("JLArrays"))
HAS_JLARRAYS && @eval using JLArrays

@testset "+ $(basename(@__FILE__))" verbose = true begin
    EP, UB = ElastoPlasm, ElastoPlasm.UnifiedBackend
    if !HAS_JLARRAYS
        @info "JLArrays is not installed: the device tests are skipped (Pkg.test() installs it)"
    else
        # JLArrays has no `synchronize` method: its kernels run when launched
        hasmethod(KernelAbstractions.synchronize, Tuple{JLBackend}) || @eval KernelAbstractions.synchronize(::JLBackend) = nothing
        UB.register!(UB.Platform(; name="JLArrays (a CPU stand-in)", backend=JLBackend(), arraytype=JLArray))
        JLArrays.allowscalar(false)
        try
            noplot = (; plot=(; status=false))
            L      = [64.1584, 64.1584/4]
            function short(c; te=0.3)                                    # the same case, its stages shortened
                s = deepcopy(c.settings)
                merge!(s.run, Dict("te" => te, "tg" => te, "tep" => te/2))
                return Case(c.layout, s)
            end
            # every array of a problem's points, mesh and basis, by path (plain values, such as SVectors, left out)
            function arrays(x, path="", out=Pair{String,Any}[])
                if isbits(x)
                    nothing
                elseif x isa AbstractArray
                    push!(out, path => x)
                elseif isstructtype(typeof(x)) && parentmodule(typeof(x)) === EP || x isa Tuple
                    foreach(f -> arrays(getfield(x, f), "$path.$f", out), fieldnames(typeof(x)))
                end
                return out
            end
            # host and device arrays agree: integers and Bools exactly, floats to round-off (relative, whole array)
            function numbers(x)
                T = eltype(x)
                T <: Number && return vec(x)
                F = fieldtype(T, 1)
                return reinterpret(F <: Number ? F : eltype(F), vec(x))   # SVector, SMatrix, stresses, strains, models
            end
            function same(x, y)
                a, b = numbers(x), numbers(y)
                eltype(a) <: AbstractFloat || return a == b
                rtol = eltype(a) == Float32 ? 1e-5 : 1e-9
                return EP.norm(a .- b) ≤ rtol * max(EP.norm(a), EP.norm(b))
            end
            configs = [
                "bsmpm, std, DP, Hencky; nonlocal, MUSL, F-bar" => (;),
                "GIMP, APIC, von Mises, hypoelastic"           => (; basis=(; which="gimpm", trsfr="apic"), material=(; plastic="VM", elastic="hypoelastic")),
                "MLS, TPIC, DP, Hencky"                        => (; basis=(; which="mlsmpm", trsfr="tpic")),
                "linear, std, von Mises; no stabilisation"     => (; basis=(; which="smpm"), material=(; plastic="VM"), nonloc=(; status=false), stab=(; musl=false, locking=false)),
                "Float32"                                      => (; dtype=(; precision="32-bit precision (or single precision)")),
            ]
            for (name, opts) ∈ configs
                @testset "$name: device run ≈ host run" begin
                    hc = short(slump(L, [16, 4]; noplot..., opts...))
                    dc = short(slump(L, [16, 4]; noplot..., opts..., backend=(; select="device")))
                    host, dev = Problem(hc, Solver(hc)), Problem(dc, Solver(dc))
                    @test EP.exec_backend(dev.solver) isa JLBackend && EP.exec_backend(host.solver) isa CPU
                    for w ∈ (elastodynamic!, elastoplastic!)
                        w(host); w(dev)
                    end
                    a, b = arrays((host.mpts, host.mesh, host.basis)), arrays((dev.mpts, dev.mesh, dev.basis))
                    @test first.(a) == first.(b)
                    @test all((((_, x), (_, y)),) -> same(x, y), zip(a, b))
                    @test all(((_, y),) -> y isa Array, b)                                  # the host problem stays on the host
                    @test dev.time.t == host.time.t
                end
            end
            @testset "3D" begin
                hc = short(slump([L; L[2]], [8, 2, 2]; noplot...); te=0.1)
                dc = short(slump([L; L[2]], [8, 2, 2]; noplot..., backend=(; select="device")); te=0.1)
                host, dev = Problem(hc, Solver(hc)), Problem(dc, Solver(dc))
                foreach(w -> (w(host); w(dev)), (elastodynamic!, elastoplastic!))
                @test all((((_, x), (_, y)),) -> same(x, y), zip(arrays((host.mpts, host.mesh)), arrays((dev.mpts, dev.mesh))))
            end
            @testset "the device copy: device arrays only, back to the same types" begin
                c = slump(L, [16, 4]; noplot..., basis=(; which="gimpm", trsfr="apic"), backend=(; select="device"))
                p = Problem(c, Solver(c))
                d = EP.to_device((p.mpts, p.mesh, p.basis), p.solver)
                @test all(((_, x),) -> x isa JLArray, arrays(d))                           # no host array left in a kernel argument
                @test typeof(EP.Adapt.adapt(Array, d)) == typeof((p.mpts, p.mesh, p.basis))
                @test EP.to_host!((p.mpts, p.mesh, p.basis), d) === (p.mpts, p.mesh, p.basis)
                @test_throws Exception d[1].x[1]                                           # element access from the host refused
            end
            @testset "dynamic relaxation and thermal refuse a device" begin
                bc = block([8, 4]; noplot..., backend=(; select="device"))
                q  = Problem(bc, Solver(bc))
                @test_throws ArgumentError elastoquasistatic!(q; lstps=1)
                tc = short(thermal(L, [8, 2]; noplot..., backend=(; select="device")))
                t  = Problem(tc, Solver(tc))
                @test_throws ArgumentError thermodynamic!(t)
            end
        finally
            JLArrays.allowscalar(true)
            filter!(p -> !(p.backend isa JLBackend), UB.device_platforms())
        end
    end
end
