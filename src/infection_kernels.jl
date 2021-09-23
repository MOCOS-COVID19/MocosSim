function enqueue_transmissions!(state::SimState, ::Val{ConstantKernelContact}, source_id::Integer, params::SimParams)
  progression = progressionof(params, source_id)

  start_time = progression.incubation_time
  end_time =  if      !ismissing(progression.mild_symptoms_time);   progression.mild_symptoms_time
              elseif  !ismissing(progression.severe_symptoms_time); progression.severe_symptoms_time
              elseif  !ismissing(progression.recovery_time);        progression.recovery_time
              else    error("no recovery nor symptoms time defined")
              end

  strain = strainof(state, source_id)

  total_infection_rate = (end_time - start_time) * params.constant_kernel_param
  total_infection_rate *= straindata(params, strain).outside_multiplier
  total_infection_rate *= spreading(params, source_id)

  num_infections = rand(state.rng, Poisson(total_infection_rate))

  if num_infections == 0
    return
  end
  @assert start_time != end_time "pathologicaly short time for infections there shouldn't be any infections but are $num_infections, progression=$progression"

  time_dist = Uniform(state.time, end_time - start_time + state.time) # in global time reference frame

  num_individuals = numindividuals(params)

  for _ in 1:num_infections
    subject_id = sample(state.rng, 1:(num_individuals-1)) # exclude the source itself
    if subject_id >= source_id # restore the indexing
      subject_id +=1
    end

    if Healthy == health(state, subject_id)
      infection_time::TimePoint = rand(state.rng, time_dist) |> TimePoint
      @assert state.time <= infection_time <= (end_time-start_time + state.time)
      push!(state.queue, Event(Val(TransmissionEvent), infection_time, subject_id, source_id, ConstantKernelContact, strain))
    end
  end
  nothing
end

function enqueue_transmissions!(state::SimState, ::Val{HouseholdContact}, source_id::Integer, params::SimParams)
  household = householdof(params, source_id)

  if 1 == length(household)
    return
  end

  progression = progressionof(params, source_id)

  start_time = progression.incubation_time
  end_time =  if      !ismissing(progression.severe_symptoms_time); progression.severe_symptoms_time
              elseif  !ismissing(progression.recovery_time);        progression.recovery_time
              else    error("no recovery nor severe symptoms time defined")
              end

  max_time = time(state) - start_time + end_time

  strain = strainof(state, source_id)
  mean_infection_time = (length(household)-1) / params.household_kernel_param
  mean_infection_time /= straindata(params, strain).household_multiplier
  time_dist = Exponential(mean_infection_time)

  for subject_id in household
    if subject_id == source_id || Healthy != health(state, subject_id)
      continue
    end

    infection_time = time(state) + rand(state.rng, time_dist)
    if infection_time > max_time
      continue
    end

    @assert time(state) <= infection_time <= (end_time - start_time + time(state))
    push!(state.queue, Event(Val(TransmissionEvent),
      infection_time,
      subject_id,
      source_id,
      HouseholdContact,
      strain)
    )
  end
  nothing
end

function enqueue_transmissions!(state::SimState, ::Val{HospitalContact}, source_id::Integer, params::SimParams)
  if nothing === params.hospital_kernel_params
    return
  end
  progression = progressionof(params, source_id)

  @assert progression.severity in SA[Severe, Critical]
  @assert !ismissing(progression.severe_symptoms_time)
  @assert !ismissing(progression.recovery_time) || !ismissing(progression.death_time)

  start_time = progression.severe_symptoms_time
  end_time = ismissing(progression.recovery_time) ? progression.death_time : progression.recovery_time

  strain = strainof(state, source_id)

  total_infection_rate = (end_time - start_time) * params.hospital_kernel_params.kernel_constant
  total_infection_rate *= straindata(params, strain).outside_multiplier

  num_infections = rand(state.rng, Poisson(total_infection_rate))

  if num_infections == 0
    return
  end
  @assert start_time != end_time "pathologicaly short time for infections there shouldn't be any infections but are $num_infections, progression=$progression"

  time_dist = Uniform(state.time, end_time - start_time + state.time) # in global time reference frame

  for _ in 1:num_infections
    subject_id = sample(state.rng, params.hospital_kernel_params.hospital_staff_ids)
    if subject_id == source_id # self infection not possible
      continue
    end

    if Healthy == health(state, subject_id)
      infection_time::TimePoint = rand(state.rng, time_dist) |> TimePoint
      @assert state.time <= infection_time <= (end_time-start_time + state.time)
      push!(state.queue, Event(Val(TransmissionEvent), infection_time, subject_id, source_id, HospitalContact))
    end
  end
end

function enqueue_transmissions!(state::SimState, ::Val{AgeCouplingContact}, source_id::Integer, params::SimParams)
  if params.age_coupling_params === nothing
    return
  end

  progression = progressionof(params, source_id)

  start_time = progression.incubation_time
  end_time =  if      !ismissing(progression.mild_symptoms_time);   progression.mild_symptoms_time
              elseif  !ismissing(progression.severe_symptoms_time); progression.severe_symptoms_time
              elseif  !ismissing(progression.recovery_time);        progression.recovery_time
              else    error("no recovery nor symptoms time defined")
              end

  strain = strainof(state, source_id)

  total_infection_rate = (end_time - start_time) * straindata(params, strain).outside_multiplier
  total_infection_rate *= sourceweightof(params.age_coupling_params, source_id)
  total_infection_rate *= spreading(params, source_id)

  num_infections = rand(state.rng, Poisson(total_infection_rate))

  if num_infections == 0
    return
  end
  @assert start_time != end_time "pathologicaly short time for infections there shouldn't be any infections but there are $num_infections, progression=$progression"

  time_dist = Uniform(state.time, end_time - start_time + state.time) # in global time reference frame

  for _ in 1:num_infections
    subject_id = sample(state.rng, params.age_coupling_params.coupling, source_id)

    if subject_id == source_id # self-infection difficult to avoid at sampling
      break
    end

    if Healthy == health(state, subject_id)
      infection_time::TimePoint = rand(state.rng, time_dist) |> TimePoint
      @assert state.time <= infection_time <= (end_time-start_time + state.time)
      push!(state.queue, Event(Val(TransmissionEvent), infection_time, subject_id, source_id, AgeCouplingContact, strain))
    end
  end
  nothing
end
