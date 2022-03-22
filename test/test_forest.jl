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
include("../src/enums.jl")
include("../src/event.jl")
include("../src/robin_forest.jl")

function trans(; src, subject, t=0.0)
    return Event(Val(TransmissionEvent), time=t, subject=subject, source=src, contact_kind=NoContact, strain=NullStrain)
end

@testset "ForestTesting" begin
    @testset "R count 1" begin
        include("../src/r_count.jl")
        f = RobinForest(10)
        push!(f, trans(src=1, subject=2, t=1.1))
        push!(f, trans(src=2, subject=3, t=2.1))
        push!(f, trans(src=2, subject=4, t=2.2))
        push!(f, trans(src=2, subject=5, t=2.2))
        push!(f, trans(src=5, subject=6, t=3.2))
        push!(f, trans(src=6, subject=7, t=3.3))
        @test 2 + 2 == 4

        println("backwardinfection $(backwardinfection(f, 2) )")
        println("backwardinfection subject $(  subject(backwardinfection(f, 2)) )")
        println("backwardinfection source $(  source(backwardinfection(f, 2)) )")

        for x in forwardinfections(f, 2)
            println("forwardinfections -> $x")
        end

        for x in f.inedges
            println("inedges -> $x")
        end

        infected_per_day = count_infected_per_time(f, 3)
        println("infected_per_day = $infected_per_day")


        r_per_time = count_r_per_time(f, 3)
        println("infected_per_day = $r_per_time")

    end

    @testset "Small Tree" begin
        f = RobinForest(10)
        push!(f, trans(src=1, subject=2))
        push!(f, trans(src=2, subject=3))
        push!(f, trans(src=2, subject=4))
        @test 2 + 2 == 4

        println("backwardinfection $(backwardinfection(f, 2) )")
        println("backwardinfection subject $(  subject(backwardinfection(f, 2)) )")
        println("backwardinfection source $(  source(backwardinfection(f, 2)) )")

        for x in forwardinfections(f, 2)
            println("forwardinfections -> $x")
        end

        for x in f.inedges
            println("inedges -> $x")
        end

    end

    @testset "ForestForward" begin

        f = RobinForest(10)

        e = Event(Val(TransmissionEvent), time=0.0, subject=2, source=1, contact_kind=NoContact, strain=NullStrain)

        push!(f, e)
        push!(f, trans(src=1, subject=2))
        @test 2 + 2 == 4

    end
end