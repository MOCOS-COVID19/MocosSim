#=
test_r_count:
- Julia version: 
- Author: jacek
- Date: 2022-03-07
=#

using Test
using DataStructures
using FixedPointNumbers
const PersonIdx=UInt32
const TimePoint = Fixed{Int32, 16}
include("../src/enums.jl")
include("../src/event.jl")
include("../src/robin_forest.jl")
include("../src/r_count.jl")

function trans(; src, subject, t=0.0)
    return Event(Val(TransmissionEvent), time=t, subject=subject, source=src, contact_kind=NoContact, strain=NullStrain)
end

function nan2zero(x)
    if isnan(x)
        return 0.0
        end
        return abs(x)
    end

    function max_abs_diff(arr1, arr2)
        maximum(map(nan2zero, arr1 .- arr2))
        end

        EPS = 0.00001


function test1()
    f = RobinForest(10)
    push!(f, trans(src=1, subject=2, t=1.1))  # day 1

    push!(f, trans(src=2, subject=3, t=2.1))  # day 2
    push!(f, trans(src=2, subject=4, t=2.2))  # day 2
    push!(f, trans(src=2, subject=5, t=2.2))  # day 2

    push!(f, trans(src=5, subject=6, t=3.2))  # day 3
    push!(f, trans(src=6, subject=7, t=3.3))  # day 3

    max_day = maximum(f.inedges) do e
        time(e) |> floor |> Int
    end

#         infected_per_day = count_infected_per_time(f, max_day)
#         println("Number of infected people per each day ($max_day days total) = $infected_per_day")
#         @test infected_per_day == [1, 3, 2]


    r_per_time = count_r_per_time(f, max_day)
    println("R value per each day ($max_day days total) = $r_per_time")
#     @test max_abs_diff(r_per_time, [1, 3, 0.33, 0.5] ) < EPS
end

test1()

#=
On Jacek machine it runs like:

Number of infected people per each day (3 days total) = [1, 3, 2]
R value per each day (3 days total) = [NaN, 0.75, 0.3333333333333333]
Test Summary: | Pass  Total
RCountTesting |    2      2

Process finished with exit code 0
=#