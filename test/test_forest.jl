#=
test_forest:
- Julia version: 
- Author: jacek
- Date: 2022-02-28
=#

using Test
using DataStructures
using FixedPointNumbers
const PersonIdx=UInt32
const TimePoint = Fixed{Int32, 16}
include("/mnt/data_sata/jacek/repo/mocos/MocosSim/src/enums.jl")
include("/mnt/data_sata/jacek/repo/mocos/MocosSim/src/event.jl")
include("/mnt/data_sata/jacek/repo/mocos/MocosSim/src/robin_forest.jl")

@testset "ForestTesting" begin
    @testset "ForestForward" begin

        f = RobinForest(10)

        @test 2 + 2 == 4

    end
end