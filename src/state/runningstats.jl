mutable struct RunningStats
  num_infected::Int
  num_detected::Int
  num_dead::Int
  daily_detections::Vector{UInt32}

  RunningStats() = new(0, 0, 0, UInt32[])
end

numinfected(stats::RunningStats) = stats.num_infected
numdetected(stats::RunningStats) = stats.num_detected
numdead(stats::RunningStats) = stats.num_dead
dailydetections(stats::RunningStats) = stats.daily_detections

function update!(stats::RunningStats, event::Event)
  ek = kind(event)
  if isdetection(ek)
    stats.num_detected += 1

    day = floor(UInt, time(event)) + 1
    while day > length(stats.daily_detections)
      push!(stats.daily_detections, 0)
    end
    stats.daily_detections[day] += 1

  elseif istransmission(ek)
    stats.num_infected += 1
  elseif isdeath(ek)
    stats.num_dead += 1
  end
  nothing
end

function reset!(stats::RunningStats)
  stats.num_infected = 0
  stats.num_detected = 0
  stats.num_infected = 0
  empty!(stats.daily_detections)
end

function sizehint!(stats::RunningStats, days::Integer)
  sizehint!(stats.daily_detections, days)
end