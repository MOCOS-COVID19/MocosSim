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

function trans(; src, subject)
    return Event(Val(TransmissionEvent), time=0.0, subject=subject, source=src, contact_kind=NoContact, strain=NullStrain)
end

@testset "ForestTesting" begin
    @testset "ForestForward" begin

        f = RobinForest(10)

        e = Event(Val(TransmissionEvent), time=0.0, subject=2, source=1, contact_kind=NoContact, strain=NullStrain)

        push!(f, e)
        push!(f, trans(src=1, subject=2))
        @test 2 + 2 == 4

    end
end