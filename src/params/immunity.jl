
const StrainInfectivityTable = SVector{NUM_STRAINS, Float64}

function make_infectivity_table(;base_multiplier::Real=1.0, british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5, omicron_multiplier::Real=1.7*1.5*2.0)::StrainInfectivityTable
  # needs validation with real data
  immunity = StrainInfectivityTable(base_multiplier, british_multiplier, delta_multiplier, omicron_multiplier)
  @assert all( 0 .<= immunity )
  immunity
end

struct ImmunizationEvents
  subjects::Vector{PersonIdx}
  times::Vector{TimePoint}
  times_lost_infection_immunity::Vector{TimePoint}
  times_lost_severe_immunity::Vector{TimePoint}
end


straininfectivity(table::StrainInfectivityTable, strain::StrainKind) = table[UInt(strain)]
immunited(immunity::ImmunityState) = immunity == Immunity


function make_immunity_table(state::AbstractSimState, level::Real)
  N = numindividuals(state)
  subjects = PersonIdx[]
  times = TimePoint[]
  times_lost_infection_immunity = TimePoint[]
  times_lost_severe_immunity = TimePoint[]
  for id in 1:floor(N*level)
    push!(subjects, id)
    push!(times, 0.0)
    push!(times_lost_infection_immunity, rand(state.rng)*90)
    push!(times_lost_severe_immunity, rand(state.rng)*210)
  end
  ImmunizationEvents(subjects, times, times_lost_infection_immunity, times_lost_severe_immunity)
end

# function immunize!(state::AbstractSimState, params::AbstractSimParams, immunization_thresholds::Vector{Int32},immunization_table::Matrix{Float32}, previously_infected::Vector{Float32})::Nothing
#   @info "Immunizing"
#   N = numindividuals(state)
#   @assert length(previously_infected) == 3
#   Samplers = AliasSampler[]
#   for gr in 1:length(immunization_thresholds)
#     no_immunity = 1-immunization_table[gr,1]
#     full_vaccinated = immunization_table[gr,1]-immunization_table[gr,2]
#     boostered = immunization_table[gr,2]
#     immunization_prob = [no_immunity*(1-previously_infected[1]), no_immunity*previously_infected[1], full_vaccinated*(1-previously_infected[2]), full_vaccinated*previously_infected[2], boostered*(1-previously_infected[3]), boostered*previously_infected[3]]
#     push!(Samplers, MocosSim.AliasSampler(UInt8,immunization_prob))
#   end
#   for id in 1:N
#     age = params.ages[id]
#     group_ids = agegroup(immunization_thresholds,age)
#     immunity_int = asample(Samplers[group_ids])
#     immunity = immunity_int |> ImmunityState
#     setimmunity!(state, id, immunity)
#   end
#   nothing
# end

function immunize!(state::SimState, immunization::ImmunizationEvents; enqueue::Bool)::Nothing
  @info "Immunizing"
  N = length(immunization.times)
  current_time = time(state)
  count = 0
  enqueued = 0
  for i in 1:N
    if immunization.times[i] <= current_time
      setimmunity!(state, immunization.subjects[i], Immunity, immunization.times_lost_infection_immunity[i], immunization.times_lost_severe_immunity[i])
      count += 1
    elseif enqueue
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.times[i],
          immunization.subjects[i],
          Immunity#,
          #immunization.times_lost_infection_immunity[i],
          #immunization.times_lost_severe_immunity[i]
        )
      )
      enqueued += 1
    end
  end

  @info "Executed $N entires: immunized $count individuals and enqueued $enqueued "

  nothing
end
