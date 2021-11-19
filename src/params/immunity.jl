const StrainImmunityTable = SMatrix{NUM_IMMUNITIES, NUM_STRAINS, Float64, NUM_IMMUNITIES*NUM_STRAINS }

function make_infectivity_table(;base_multiplier::Real=1.0, british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5)::StrainImmunityTable
  # needs validation with real data
  mat = @SMatrix [
    1.00    1.70    2.55;
    0.01    0.01    0.10;
    0.10    0.10    0.50;
    0.03    0.10    0.30;
  ]

  @assert all( mat .>= 0)
  @assert maximum(mat, dims=1) |> vec == mat[1,:]
  mat
end

struct Immunization
  times::Vector{TimePoint}
  subjects::Vector{PersonIdx}
  immunity_kinds::Vector{ImmunityState}
end

rawinfectivity(table::StrainImmunityTable, strain::StrainKind) = table[UInt(1), UInt(strain)]
condinfectivity(table::StrainImmunityTable, immunity::ImmunityState, strain::StrainKind) = table[UInt(immunity), UInt(strain)] / table[UInt(1), UInt(strain)]

function immunize!(state::SimState, immunization::Immunization; enqueue::Bool)::Nothing
  N = length(immunization.times)
  @assert N == length(immunization.subjects) == length(immunization.immunity_kinds)
  @assert issorted(immunization.times)

  current_time = time(state)
  i = 1

  while i <= N && current_time <= immunization.times[i]
    setimmunity!(state, immunization.subjects[i], immunization.immunity_kinds[i])
  end

  if !enqueue
    return nothing
  end

  while i <= N
    push!(state.queue,
      Event(
        Val(ImmunizationEvent),
        immunization.times[i],
        immunization.subjects[i],
        immunization.immunity_kinds[i],
      )
    )
  end

  nothing
end