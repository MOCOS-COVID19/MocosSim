using MocosSim
using Test
using Random

tests = [
  "population_grouping"
  "matrix_alias_sampler"
]

if length(ARGS) > 0
    tests = ARGS
end

@testset "MocosSim" begin

for t in tests
    fp = joinpath(dirname(@__FILE__), "test_$t.jl")
    println("$fp ...")
    include(fp)
end

end # @testset