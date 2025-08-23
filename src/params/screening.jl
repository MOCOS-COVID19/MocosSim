Base.@kwdef struct ScreeningParams
  start_time::Float64 = 0.0
  precision::Float64 = 0.8
  period::Float64 = 7.0
  lower_bound_age::Int64 = 8
  upper_bound_age::Int64 = 16
  test_times::Vector{TimePoint} = TimePoint[]
  adherence_pmf::Vector{Float64} = Float64[]  # length should be 86
end

function screening!(state::AbstractSimState, params::AbstractSimParams, event::Event)
  if params.screening_params != nothing
    # --- CACHE: stores school_id → Bool (true = participates, false = skips)
    school_participation = Dict{Int, Bool}()
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
      # --- Determine school participation only ONCE ---
      school_id = school(params, id)
      if !haskey(school_participation, school_id)
          school_adherence_probability = school_adherence_probability(params, id)
          school_participation[school_id] = rand(state.rng) <= school_adherence_probability
      end
      # --- Skip the entire school if it opted out ---
      if school_participation[school_id] == false
          continue
      end
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

function add_screening!(state::AbstractSimState, params::AbstractSimParams)
  if !isempty(params.screening_params.test_times)
    for screening_time in params.screening_params.test_times
      event = Event(Val(ScreeningEvent), screening_time)
      push!(state.queue, event)
    end
  end
end