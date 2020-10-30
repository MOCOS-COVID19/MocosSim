struct PhoneTracingParams
  isusingapp::BitVector
  detection_prob::Float64
  detection_delay::Float64
end

PhoneTracingParams(N::Integer, detection_prob::Real=1.0, detection_delay::Real=0.5) =
  PhoneTracingParams(BitVector(undef, N), detection_prob, detection_delay)

function resample!(rng::AbstractRNG, params::PhoneTracingParams, prob::Real)
  for i in 1:length(params.isusingapp)
    params.isusingapp[i] = rand(rng) < prob
  end
  params
end

function resamplebyhouseholds!(rng::AbstractRNG, params::PhoneTracingParams, prob::Real, household_ptrs::Vector{Tuple{Ti,Ti}} where Ti<:Integer)
  num_individuals = length(household_ptrs)
  num_users = rand(rng, Binomial(num_individuals, prob)) # makes the number same as in the independent case

  randomized_households = shuffle!(rng, unique(household_ptrs)) # not optimal but gets the job done

  fill!(params.isusingapp, false)
  num_assigned = 0
  for (household_start, household_end) in randomized_households
    if num_assigned >= num_users
      break
    end

    household = UnitRange(household_start, household_end)

    assignements_in_household = min(length(household), num_users-num_assigned)
    for member in 1:assignements_in_household
      params.isusingapp[household[member]] = true
    end
    num_assigned += assignements_in_household
  end
end

function PhoneTracingParams(rng::AbstractRNG,
                             N::Integer,
                             usage::Real,
                             detection_delay::Real=0.25,
                             detection_prob::Real=1.0,
                             household_ptrs::Union{Nothing,Vector{Tuple{Ti,Ti}} where Ti <: Integer} = nothing)
  params = PhoneTracingParams(N, detection_prob, detection_delay)
  if household_ptrs === nothing
    resample!(rng, params, usage)
  else
    resamplebyhouseholds!(rng, params, usage, household_ptrs)
  end
  params
end

uses_phone_tracing(params::PhoneTracingParams, person_id::Integer) = params.isusingapp[person_id]

function saveparams(dict, p::PhoneTracingParams, prefix::AbstractString="")
  dict[prefix*"isusingapp"] = p.isusingapp
  dict[prefix*"detection_prob"] = p.detection_prob
  dict[prefix*"detection_delay"] = p.detection_delay
end
