using Test,Plots,LaTeXStrings,Revise


using ElastoPlasm

@testset "ElastoPlasm.jl" verbose = true begin
    @testset "└ shpTest.jl" verbose = true begin
        shpTest(;ghost=true)
    end
    @testset "└ slumpTest.jl" verbose = true begin
        slumpTest()
    end
end