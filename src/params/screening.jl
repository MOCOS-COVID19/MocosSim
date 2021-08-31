Base.@kwdef struct ScreeningParams <: ScreeningParam
  start_time::Float64 = 0.0
  precision::Float64 = 0.8
  period::Float64 = 7.0
  lower_bound_age::Int64 = 8
  upper_bound_age::Int64 = 16
end

function screening!(state::SimState, params::SimParams, event::Event)
  if params.screening_params != nothing
    for id in 1:numindividuals(state)
      health = MocosSim.health(state, id)
      if health == Healthy || health == Recovered || health == Incubating
        continue
      end

      age = params.ages[id]
      if age < params.screening_params.lower_bound_age || age > params.screening_params.upper_bound_age
        continue
      end

      if rand(state.rng) >= params.screening_params.precision
        continue
      end
      push!(
        state.queue, 
        Event(
          Val(DetectionEvent),
          time(event),
          id,
          OutsideQuarantineDetction),
          immediate=true)
    end
  end
end

function add_screening!(state::SimState, params::SimParams,time_limit::TimePoint=typemax(TimePoint))
  for screening_time in 0.0:params.screening_params.period:time_limit
    event = Event(Val(ScreenigEvent), screening_time)
    push!(state.queue, event)
  end
end
