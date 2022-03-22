#=
r_count:
- Julia version: 
- Author: jacek
- Date: 2022-03-07
=#
using DataStructures: Deque
using Debugger
include("../src/default_vector.jl")

function count_infected_per_time(forest::RobinForest, max_days::Int, time_unit=1.0)
    infected_per_time = Array{Int}(undef, max_days)
    fill!(infected_per_time, 0)
    N = size(forest.inedges)
    for e in forest.inedges
        if kind(e) == TransmissionEvent
            day_number = (time(e) / time_unit) |> floor |> Int
            infected_per_time[day_number] += 1
            end
        end
        infected_per_time
    end

function count_r_per_time(forest::RobinForest, max_days::Int, time_unit=1.0)
    r_per_time = DefaultVector(-1.0, Float64)
    infected_per_time = DefaultVector(0, Int)
    edge_springsoff_per_time = DefaultVector(0, Int)

    for e in forest.inedges
        println(e)
        end
        println()
        println("outdegrees")
    for e in forest.outdegrees
        println(e)
        end
        println()
        println("outedgedict")
    for e in forest.outedgedict
        println(e)
        end

    for p in forwardinfections(forest, 2)
        println("2 -> $(subject(p)) $p")
        end


    init_node = 1
    d = Deque{Int64}()
    pushfirst!(d, 1)
    pushfirst!(d, 2)

    println(pop!(d))
    println(pop!(d))

    person_infecting = zeros(Int, length(forest.inedges))

    person_infected_at = Vector{Int}(undef, length(forest.inedges))
    person_infected_at[1] = -1
    person_infecting[1] = 1
    pushfirst!(d, 1)
    @bp

    while !isempty(d)
        n = pop!(d)
        for e::Event in forwardinfections(forest, n)
            n2 = subject(e)
            pushfirst!(d, n2)
            t = (time(e) / time_unit) |> floor |> Int
            println("n = $n, n2 = $n2, t = $t")
            infected_per_time[t] += 1
            person_infected_at[n2] = t
            person_infecting[n] = t
#             if person_infected_at[n] + 1 == person_infected_at[n2]
#                 edge_springsoff_per_time[t] += 1
                edge_springsoff_per_time[t] += 1
#             end
            println("person_infected_at[$t] = $(person_infected_at[t])")
        end
    end

    person_infecting_at_time = zeros(Int, length(forest.inedges))
    for i in eachindex(person_infecting)
        if person_infecting[i] > 0
                person_infecting_at_time[person_infecting[i]] += 1
            end
    end


#     return convert(Vector, r_per_time)
    println("edge_springsoff_per_time = $(convert(Vector, edge_springsoff_per_time))")
    println("infected_per_time = $(convert(Vector, infected_per_time))")
    println("person_infecting = $(convert(Vector, person_infecting))")
    println("person_infecting_at_time = $(convert(Vector, person_infecting_at_time))")
    return [(springs_off / per_day) for (springs_off, per_day) in zip(convert(Vector, infected_per_time), convert(Vector, person_infecting_at_time))]
    end