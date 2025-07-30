push!(LOAD_PATH, "../src")

using Test,JLD2,ProgressMeter,Suppressor,Plots,LaTeXStrings,ElastoPlasm
using BenchmarkTools, KernelAbstractions
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync

function runtests()
    global ROOT    = dirname(@__FILE__)
    global DATASET = joinpath(ROOT, "dataset")
    if !isdir(DATASET)
        mkpath(DATASET)
    end

    testdir   = joinpath(@__DIR__,"testset")
    istest(f) = endswith(f, ".jl") && startswith(f, "test_")
    testfiles = sort(filter(istest, readdir(testdir)))
    nfail = 0
    testfiles = ["test_performance.jl"]
    @testset "ElastoPlasm.jl tested:" verbose = true begin
        for f ∈ testfiles
            try
                printstyled("Running: $(f)\n"; bold=true, color=:white)
                include(joinpath(testdir, f))
            catch ex
                nfail += 1
            end
            println("")
        end
    end
    return nfail
end
exit(runtests())

#=
To reproduce CI.yaml locally, run the following from the same repository state on julia version 1.10.10:
    julia> import Pkg

    julia> Pkg.test(;coverage=true, julia_args=["--check-bounds=yes", "--compiled-modules=yes"], force_latest_compatible_version=false, allow_reresolve=true)
=#