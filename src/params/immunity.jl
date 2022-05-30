
const StrainInfectivityTable = SVector{NUM_STRAINS, Float64}

function make_infectivity_table(;base_multiplier::Real=1.0, british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5, omicron_multiplier::Real=1.7*1.5*2.0)::StrainInfectivityTable
  # needs validation with real data
  immunity = StrainInfectivityTable(base_multiplier, british_multiplier, delta_multiplier, omicron_multiplier)
  @assert all( 0 .<= immunity )
  immunity
end

struct ImmunizationEvents
  subjects::Vector{PersonIdx}
  immunity_uptake_times::Vector{TimePoint}
  lost_immunity_times::Vector{TimePoint}
  immunity_types::Vector{ImmunityState}
end

function make_immunity_table(state::AbstractSimState, level::Real)
  N = numindividuals(state)
  subjects = PersonIdx[]
  immunity_uptake_times = TimePoint[]
  lost_immunity_times = TimePoint[]
  immunity_types = ImmunityState[]
  for id in 1:floor(N*level)
    push!(subjects, id)
    push!(immunity_uptake_times, rand(state.rng)*30)
    push!(lost_immunity_times, rand(state.rng)*90)
    push!(immunity_types, against_infection)
  end
  ImmunizationEvents(subjects, immunity_uptake_times, lost_immunity_times, immunity_types)
end

straininfectivity(table::StrainInfectivityTable, strain::StrainKind) = table[UInt(strain)]
immunited(immunity::ImmunityState) = immunity != NoImmunity


function immunize!(state::SimState, immunization::ImmunizationEvents)::Nothing
  @info "Immunizing"
  N = length(immunization.immunity_uptake_times)
  current_time = time(state)
  count = 0
  enqueued = 0
  for i in 1:N
    if immunization.immunity_uptake_times[i] <= current_time
      setimmunity!(state, immunization.subjects[i], immunization.immunity_types[i])
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.lost_immunity_times[i],
          immunization.subjects[i],
          NoImmunity
        )
      )
      count += 1
    else
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.immunity_uptake_times[i],
          immunization.subjects[i],
          immunization.immunity_types[i]
        )
      )
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.lost_immunity_times[i],
          immunization.subjects[i],
          NoImmunity
        )
      )
      enqueued += 1
    end
  end

  @info "Executed $N entires: immunized $count individuals and enqueued $enqueued "


  nothing
end
