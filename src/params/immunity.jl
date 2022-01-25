
const StrainInfectivityTable = SVector{NUM_STRAINS, Float64}
const StrainSusceptibilityTable = SMatrix{NUM_IMMUNITIES, NUM_STRAINS, Float64, NUM_IMMUNITIES*NUM_STRAINS }

function make_infectivity_table(;base_multiplier::Real=1.0, british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5, omicron_multiplier::Real=1.7*1.5*2.0)::StrainInfectivityTable
  # needs validation with real data
  immunity = StrainInfectivityTable(base_multiplier, british_multiplier, delta_multiplier, omicron_multiplier)
  @assert all( 0 .<= immunity )
  immunity
end

function make_susceptibility_table()
  #each column is distinct StrainKind
  #each row is distinct ImmunityState

  table = @SMatrix [
    #ChineseStrain  BritishStrain DeltaStrain OmicronStrain
    1.00            1.00          1.00        1.00;         # NoImmunity
    0.01            0.01          0.10        0.60;         # NaturalImmunity
    0.10            0.10          0.50        0.75;         # VecVacImmunity
    0.03            0.10          0.30        0.75;         # MRNAVacImmunity
  ]

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
