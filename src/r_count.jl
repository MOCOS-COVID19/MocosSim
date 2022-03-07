#=
r_count:
- Julia version: 
- Author: jacek
- Date: 2022-03-07
=#

function count_infected_per_day(forest::RobinForest, max_days::Int)
    infected_per_day = Array{Int}(undef, max_days)
    fill!(infected_per_day, 0)
    N = size(forest.inedges)
    for e in forest.inedges
        if kind(e) == TransmissionEvent
            day_number = time(e) |> floor |> Int
            infected_per_day[day_number] += 1
            end
        end
        infected_per_day
    end