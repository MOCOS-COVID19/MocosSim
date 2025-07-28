Base.@kwdef struct ScreeningParams
  start_time::Float64 = 0.0
  precision::Float64 = 0.8
  period::Float64 = 7.0
  lower_bound_age::Int64 = 8
  upper_bound_age::Int64 = 16
  interval_periods::Vector{Float64} = Float64[]
  interval_times::Vector{TimePoint} = TimePoint[]
end

function screening!(state::AbstractSimState, params::AbstractSimParams, event::Event)
  if params.screening_params != nothing
    for id in 1:numindividuals(state)
      health = MocosSim.health(state, id)
      # quick speedup as test says No in these scenarios
      if health == Healthy || health == Recovered || health == Incubating
        continue
      end
      # age condition:
      age = params.ages[id]
      if age < params.screening_params.lower_bound_age || age > params.screening_params.upper_bound_age
        continue
      end
      # test precision condition:
      if rand(state.rng) >= params.screening_params.precision
        continue
      end
      screening_freedom = freedom(state, id)
      # school eligibility check conditions:
      if ((HomeTreatment == screening_freedom) || (HomeQuarantine == screening_freedom) || (Hospitalized == screening_freedom))
          continue
      end
      if
      push!(
        state.queue, 
        Event(
          Val(DetectionEvent),
          time(event),
          id,
          OutsideQuarantineScreeningDetection),
          immediate=true)
    end
  end
end

function add_screening!(state::AbstractSimState, params::AbstractSimParams, time_limit::TimePoint=typemax(TimePoint))
  if length(params.screening_params.interval_times) == 0
    for screening_time in params.screening_params.start_time:params.screening_params.period:time_limit
      event = Event(Val(ScreeningEvent), screening_time)
      push!(state.queue, event)
    end
  else
    t = copy(params.screening_params.interval_times)
    periods = params.screening_params.interval_periods
    @assert length(t) == length(periods)
    append!(t, time_limit)
    for index in 1 : length(t) - 1
      last_elem = missing
      for ti in t[index] : periods[index] : t[index + 1]
        last_elem = ti
        event = Event(Val(ScreeningEvent), ti)
        push!(state.queue, event)
      end
      if !ismissing(last_elem)
        if index < length(t) - 1
          t[index + 1] = last_elem + periods[index + 1]
        end
      end
    end
  end
end