# Validate ElastoPlasm's explicit solver on an NVIDIA GPU (CUDA), against the host. From the repository root, on the
# GPU machine (cloned with --recursive):
#
#     julia -t 1 test/cuda/validate_cuda.jl
#
# It sets up its own environment in test/cuda (ElastoPlasm, EasyDraw and UnifiedBackend from this checkout, plus
# CUDA), then: reports the device; compiles and runs a few steps of every configuration, collecting all errors;
# compares device and host runs; times a step on both; checks that a simulation saved there loads without CUDA. It
# prints a summary and writes dump/cuda_validation/report.toml (send that file back). To try the script without a GPU,
# with JLArrays as the device (small sizes only): EP_DEVICE=jlarrays julia -t 1 test/cuda/validate_cuda.jl
using Pkg
const ROOT = dirname(dirname(@__DIR__))
const MODE = get(ENV, "EP_DEVICE", "cuda")
Pkg.activate(@__DIR__; io=devnull)
Pkg.develop([PackageSpec(path=joinpath(ROOT, "lib", "EasyDraw")), PackageSpec(path=joinpath(ROOT, "lib", "UnifiedBackend")),
             PackageSpec(path=ROOT)]; io=devnull)
Pkg.add(MODE == "cuda" ? "CUDA" : "JLArrays"; io=devnull)
using TOML, LinearAlgebra
MODE == "cuda" ? (@eval using CUDA) : (@eval using JLArrays)
using ElastoPlasm
const EP, UB, KA = ElastoPlasm, ElastoPlasm.UnifiedBackend, ElastoPlasm.KernelAbstractions
report = Dict{String,Any}("mode" => MODE, "julia" => string(VERSION), "threads" => Threads.nthreads())

# ---- 1. the device
if MODE == "cuda"
    report["device"] = Dict("name" => CUDA.name(CUDA.device()), "capability" => string(CUDA.capability(CUDA.device())),
                            "cuda_driver" => string(CUDA.driver_version()), "cuda_runtime" => string(CUDA.runtime_version()),
                            "CUDA.jl" => string(pkgversion(CUDA)))
    CUDA.capability(CUDA.device()) ≥ v"6.0" || @warn "compute capability below 6.0: no native Float64 atomic add"
else
    @eval KA.synchronize(::JLBackend) = nothing                                   # JLArrays has none
    UB.register!(UB.Platform(; name="JLArrays", backend=JLBackend(), arraytype=JLArray))
    report["device"] = Dict("name" => "JLArrays (CPU stand-in)")
end
isempty(UB.device_platforms()) && error("no GPU registered by UnifiedBackend: is CUDA functional? (CUDA.functional())")
const DEVICE = (; backend=(; select="device"))
noplot = (; plot=(; status=false))
L2 = [64.1584, 64.1584/4]
with_run(c; te, tep=te/2) = Case(c.layout, Settings(; zones=c.settings.zones, defaults=c.settings.defaults, sides=c.settings.sides,
                                                       options=c.settings.options, run=merge(c.settings.run, Dict("te" => te, "tg" => te, "tep" => tep))))
println("device: ", report["device"]["name"], ", backend ", EP.exec_backend(get_solver(; DEVICE...)))

# ---- 2. every configuration compiles and runs a few steps
configs = [
    "defaults: bsmpm, std, DP, Hencky, nonlocal, MUSL, F-bar" => (;),
    "GIMP" => (; basis=(; which="gimpm")), "uGIMP" => (; basis=(; which="gimpm", domain="undeformed")),
    "linear (smpm)" => (; basis=(; which="smpm")), "MLS" => (; basis=(; which="mlsmpm")),
    "TPIC" => (; basis=(; trsfr="tpic")), "APIC" => (; basis=(; trsfr="apic")), "GIMP, APIC" => (; basis=(; which="gimpm", trsfr="apic")),
    "von Mises" => (; material=(; plastic="VM")), "hypoelastic" => (; material=(; elastic="hypoelastic")),
    "improved Hencky" => (; material=(; elastic="improved_hencky")),
    "no nonlocal, MUSL, F-bar" => (; nonloc=(; status=false), stab=(; musl=false, locking=false)),
    "Float32" => (; dtype=(; precision="32-bit precision (or single precision)")),
    "Float32, von Mises" => (; dtype=(; precision="32-bit precision (or single precision)"), material=(; plastic="VM")),
]
compile = Dict{String,Any}()
for (dims, nel) ∈ (("2D", [16, 4]), ("3D", [8, 2, 2])), (name, opts) ∈ configs
    dims == "3D" && name ∉ ("defaults: bsmpm, std, DP, Hencky, nonlocal, MUSL, F-bar", "GIMP", "MLS", "APIC", "von Mises",
                            "hypoelastic", "Float32") && continue
    key = "$dims, $name"
    t = @elapsed compile[key] = try
        L = dims == "2D" ? L2 : [L2; L2[2]]
        p = Problem(with_run(slump_case(L, nel; noplot..., opts..., DEVICE...); te=0.02))
        elastodynamic!(p); elastoplastic!(p)
        all(v -> all(isfinite, v), p.mpts.x) && all(isfinite, p.mpts.J) ? "ok" : "ok (ran; non-finite values)"
    catch e
        first(sprint(showerror, e), 2000)
    end
    println(rpad(key, 64), startswith(compile[key], "ok") ? compile[key] : "ERROR", "  ($(round(t; digits=1)) s)")
end
report["compile"] = compile

# ---- 3. device vs host: the same runs, relative differences (whole arrays) and totals. Plasticity amplifies
# round-off where a point sits at the yield threshold, so in a plastic run single fields (the plastic strain above all;
# percents in Float32) can differ much more than the totals (run-out, largest plastic strain, kinetic energy).
function fields(p)
    m = p.mpts
    num(v) = Float64.(collect(reinterpret(eltype(eltype(v)), v)))
    return Dict("x" => num(m.x), "v" => num(m.s.v), "J" => Float64.(m.J), "p" => Float64.(getfield.(m.s.σᵢ, :p)),
                "dev" => num(getfield.(m.s.σᵢ, :dev)), "epII" => num(m.s.ϵpII))
end
totals(p) = Dict("runout" => Float64(maximum(getindex.(p.mpts.x, 1))), "epII" => Float64(maximum(first.(p.mpts.s.ϵpII))),
                 "kinetic" => Float64(sum(p.mpts.s.m .* sum.(abs2, p.mpts.s.v)) / 2))
accuracy = Dict{String,Any}()
for (name, opts, workflows, te) ∈ (("elastic, Float64", (;), (elastodynamic!,), 0.5),
                                   ("plastic, DP, Float64", (;), (elastodynamic!, elastoplastic!), 1.0),
                                   ("plastic, von Mises, Float64", (; material=(; plastic="VM")), (elastodynamic!, elastoplastic!), 1.0),
                                   ("plastic, DP, Float32", (; dtype=(; precision="32-bit precision (or single precision)")), (elastodynamic!, elastoplastic!), 1.0))
    try
        host = Problem(with_run(slump_case(L2, [40, 10]; noplot..., opts...); te))
        dev  = Problem(with_run(slump_case(L2, [40, 10]; noplot..., opts..., DEVICE...); te))
        foreach(w -> (w(host); w(dev)), workflows)
        a, b = fields(host), fields(dev)
        rel = Dict(k => norm(a[k] .- b[k]) / max(norm(a[k]), floatmin()) for k ∈ keys(a))
        th, td = totals(host), totals(dev)
        accuracy[name] = Dict("relative" => rel, "totals_host" => th, "totals_device" => td,
                              "totals_relative" => Dict(k => abs(th[k] - td[k]) / max(abs(th[k]), floatmin()) for k ∈ keys(th)))
        println(rpad(name, 32), "largest relative difference ", maximum(values(rel)), "; run-out host ", th["runout"], ", device ", td["runout"])
    catch e
        accuracy[name] = first(sprint(showerror, e), 2000); println(rpad(name, 32), "ERROR")
    end
end
report["accuracy"] = accuracy

# ---- 4. time per elastoplastic step, host (this process's threads) and device, after a warm-up
function step_time(p; steps=20, warm=3)
    m, g_, b, t, s = p.mpts, p.mesh, p.basis, p.time, p.solver
    dm, dg, db = EP.to_device((m, g_, b), s)
    function step!()
        g, dt = EP.get_spacetime(dm, dg, t, t.te)
        EP.ignite(dm, dg, db, s); EP.mapsto(dm, dg, db, g, dt, s); EP.elastoplast(dm, dg, db, dt, s)
        t.t[1] += dt
    end
    foreach(_ -> step!(), 1:warm)
    return (@elapsed foreach(_ -> step!(), 1:steps)) / steps
end
sizes = MODE == "cuda" ? [[80, 20], [160, 40], [320, 80], [40, 10, 10]] : [[16, 4]]
timing = Dict{String,Any}()
for nel ∈ sizes, (prec, dt) ∈ (("Float64", (;)), ("Float32", (; dtype=(; precision="32-bit precision (or single precision)"))))
    key = "$(join(nel, "x")), $prec"
    try
        L = length(nel) == 2 ? L2 : [L2; L2[2]]
        h = step_time(Problem(slump_case(L, nel; noplot..., dt...)))
        d = step_time(Problem(slump_case(L, nel; noplot..., dt..., DEVICE...)))
        n = Problem(slump_case(L, nel; noplot..., dt...)).mpts.nmp
        timing[key] = Dict("particles" => Int(n), "host_s_per_step" => h, "device_s_per_step" => d, "speedup" => h / d)
        println(rpad(key, 24), "$n points: host $(round(1e3h; digits=2)) ms, device $(round(1e3d; digits=2)) ms per step, ×$(round(h/d; digits=1))")
    catch e
        timing[key] = first(sprint(showerror, e), 2000); println(rpad(key, 24), "ERROR")
    end
end
report["timing"] = timing

# ---- 5. a simulation saved with a device solver loads where CUDA is not loaded (ElastoPlasm's own environment)
file = joinpath(mktempdir(), "device_simulation.jld2")
save_simulation(Simulation(Problem(slump_case(L2, [16, 4]; noplot..., DEVICE...))), file)
cmd = `$(Base.julia_cmd()) --project=$ROOT -e 'using Pkg; Pkg.instantiate(; io=devnull); using ElastoPlasm; s = load_simulation(ARGS[1]); print(nameof(typeof(ElastoPlasm.exec_backend(s.problem.solver))))' $file`
report["roundtrip"] = try
    "loaded, on " * last(split(readchomp(cmd), "\n"))   # its last line (ElastoPlasm prints a banner first)
catch e
    first(sprint(showerror, e), 2000)
end
println("saved with a device solver, loaded without CUDA: ", report["roundtrip"])

out = joinpath(ROOT, "dump", "cuda_validation", "report.toml")
mkpath(dirname(out))
open(io -> TOML.print(io, report; sorted=true), out, "w")
println("report: ", out)
