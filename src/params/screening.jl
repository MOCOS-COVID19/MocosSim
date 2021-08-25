struct ScreeningParams
    begin_day::Float64
    precision::Float64
    frequency::Float64
  end


function screen_children!(rng::AbstractRNG, state::SimState, params::SimParams)
    for screening_day in params.screening_params.begin_day:params.screening_params.frequency:365
        for id in 1:numindividuals(state)
            health = MocosSim.health(state, id)
            if health == Healthy || health == Recovered || health == Incubating
              continue
            end

            age = params.ages[id]
            if age < 6 || age > 16
              continue
            end

            if rand(rng) >= params.screening_params.precision
              continue
            end
            push!(
              state.queue, 
              Event(
                Val(DetectionEvent),
                time(state),
                id,
                OutsideQuarantineDetction),
                immediate=true)

        end
    end
end