push!(LOAD_PATH, "../src")

import ElastoPlasm

function runtests()
    exename   = joinpath(Sys.BINDIR, Base.julia_exename())
    testdir   = joinpath(pwd(),"testset")
    istest(f) = endswith(f, ".jl") && startswith(f, "test_")
    testfiles = sort(filter(istest, readdir(testdir)))

    nfail = 0
    printstyled("\nTesting ElastoPlasm.jl\n"; bold=true, color=:white)
    for f in testfiles
        println("")
        try
            run(`$exename -O3 --startup-file=no --check-bounds=no $(joinpath(testdir, f))`)
        catch ex
            nfail += 1
        end
    end
    return nfail
end
exit(runtests())


#=
│ To reproduce CI.yaml locally, run the following from the same repository state on julia version 1.10.10:
│
│ `import Pkg; Pkg.test(;coverage=true, julia_args=["--check-bounds=yes", "--compiled-modules=yes"], force_latest_compatible_version=false, allow_reresolve=true)`
=#