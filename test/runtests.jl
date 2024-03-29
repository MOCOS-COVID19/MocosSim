using MocosSim
using Test
using Random

import MocosSim: time, numdead, numdetected

tests = [
  "age_coupling",
  "matrix_alias_sampler",
  "population_grouping",
  "infection_modulations",
  "household_grouping",
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