
const StrainInfectivityTable = SVector{NUM_STRAINS, Float64}
const StrainSusceptibilityTable = SMatrix{NUM_IMMUNITIES, NUM_STRAINS, Float64, NUM_IMMUNITIES*NUM_STRAINS }

function make_infectivity_table(;base_multiplier::Real=1.0, british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5, omicron_multiplier::Real=1.7*1.5*2.0)::StrainInfectivityTable
  # needs validation with real data
  immunity = StrainInfectivityTable(base_multiplier, british_multiplier, delta_multiplier, omicron_multiplier)
  @assert all( 0 .<= immunity )
  immunity
end

function make_susceptibility_table(;delta_susceptibility::Vector{T}, omicron_susceptibility::Vector{T}) where T <: Real
  #each column is distinct StrainKind
  #each row is distinct ImmunityState

  table = @SMatrix [
    #ChineseStrain  BritishStrain DeltaStrain
    1.00   1.00;
    0.01   0.01;
    0.10   0.10;
    0.10   0.10;
    0.03   0.10;
    0.03   0.10;
  ]
  table = hcat(table, delta_susceptibility, omicron_susceptibility)
  @assert all( 0 .<= table .<= 1)
  table
end

struct ImmunizationOrder
  times::Vector{TimePoint}
  subjects::Vector{PersonIdx}
  immunity_kinds::Vector{ImmunityState}
end

straininfectivity(table::StrainInfectivityTable, strain::StrainKind) = table[UInt(strain)]
susceptibility(table::StrainSusceptibilityTable, immunity::ImmunityState, strain::StrainKind) = table[UInt(immunity), UInt(strain)]
immunited(immunity::ImmunityState) = immunity != NoImmunity
vaccinated(immunity::ImmunityState) = immunity == VacImmunity || immunity == VacNaturalImmunity || immunity == BoostImmunity || immunity == BoostNaturalImmunity
boostered(immunity::ImmunityState) = immunity == BoostImmunity || immunity == BoostNaturalImmunity
previnfected(immunity::ImmunityState) = immunity == NaturalImmunity || immunity == VacNaturalImmunity || BoostNaturalImmunity


function immunize!(state::AbstractSimState, params::AbstractSimParams, immunization_thresholds::Vector{Int32},immunization_table::Matrix{Float32}, previously_infected::Vector{Float32})::Nothing
  @info "Immunizing"
  N = numindividuals(state)
  @assert length(previously_infected) == 3
  Samplers = AliasSampler[]
  for gr in 1:length(immunization_thresholds)
    no_immunity = 1-immunization_table[gr,1]
    full_vaccinated = immunization_table[gr,1]-immunization_table[gr,2]
    boostered = immunization_table[gr,2]
    immunization_prob = [no_immunity*(1-previously_infected[1]), no_immunity*previously_infected[1], full_vaccinated*(1-previously_infected[2]), full_vaccinated*previously_infected[2], boostered*(1-previously_infected[3]), boostered*previously_infected[3]]
    push!(Samplers, MocosSim.AliasSampler(UInt8,immunization_prob))
  end
  for id in 1:N
    age = params.ages[id]
    group_ids = agegroup(immunization_thresholds,age)
    immunity_int = asample(Samplers[group_ids])
    immunity = immunity_int |> ImmunityState
    setimmunity!(state, id, immunity)
  end
  nothing
end

function immunize!(state::SimState, immunization::ImmunizationOrder; enqueue::Bool)::Nothing
  @info "Immunizing"
  N = length(immunization.times)
  @assert N == length(immunization.subjects) == length(immunization.immunity_kinds)
  @assert issorted(immunization.times)

  current_time = time(state)
  count = 0
  enqueued = 0
  for i in 1:N
    if immunization.times[i] <= current_time
      setimmunity!(state, immunization.subjects[i], immunization.immunity_kinds[i])
      count += 1
    elseif enqueue
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.times[i],
          immunization.subjects[i],
          immunization.immunity_kinds[i]
        )
      )
      enqueued += 1
    end
  end

  @info "Executed $N entires: immunized $count individuals and enqueued $enqueued "


  nothing
end
