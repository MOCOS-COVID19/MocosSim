#=
r_count:
- Julia version: 
- Author: jacek
- Date: 2022-03-07
=#

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
    r_per_time = Array{Float64}(undef, max_days)
    fill!(r_per_time, 0)

    infected_per_time = count_infected_per_time(forest, max_days, time_unit)

    r_per_time[1] = NaN64
    total_infected = infected_per_time[1]
    for i in 2:max_days
        total_infected += infected_per_time[i]
        r_per_time[i] = infected_per_time[i] / total_infected
        end
        r_per_time
    end