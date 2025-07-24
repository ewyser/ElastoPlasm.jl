push!(LOAD_PATH, "../src") # FIXME: to be removed everywhere?

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