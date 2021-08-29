using CodecZlib
using DataFrames
using FileIO
using JLD2
using MocosSim
using Test
using Random

tests = [
  "age_coupling",
  "matrix_alias_sampler",
  "population_grouping",
  "modulations"
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