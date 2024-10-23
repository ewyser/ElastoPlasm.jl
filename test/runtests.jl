using Test,Plots,LaTeXStrings,Revise


using ElastoPlasm

@testset "ElastoPlasm.jl" verbose = true begin
    @testset "└ shpTest.jl" verbose = true begin
        shpTest(;ghost=true)
    end
end