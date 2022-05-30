using Random
using Distributions
#using FunctionWrappers

const Age=UInt8
const ImmunityRand = Normed{UInt32, 32}

abstract type AbstractSimParams end

abstract type InfectionModulation end

abstract type AbstractOutsideCases end

include("params/age_coupling.jl")
include("params/households.jl")
include("params/friendship.jl")
include("params/immunity.jl")
include("params/modulations.jl")
include("params/progression.jl")
include("params/hospital.jl")
include("params/phonetracing.jl")
include("params/spreading.jl")
include("params/outside_cases.jl")
include("params/screening.jl")

struct SimParams <: AbstractSimParams
  household_ptrs::Vector{Tuple{PersonIdx,PersonIdx}}  # (i1,i2) where i1 and i2 are the indices of first and last member of the household

  ages::Vector{Age}
  genders::BitVector
  immunity_rand::Vector{ImmunityRand}

  progression_params::ProgressionParams
  strain_infectivity_table::StrainInfectivityTable

  constant_kernel_param::Float32
  household_kernel_param::Float32
  age_coupling_params::Union{Nothing, AgeCouplingParams} # nothing if kernel not active
  hospital_kernel_params::Union{Nothing, HospitalInfectionParams}  # nothing if hospital kernel not active

  hospital_detections::Bool
  mild_detection_prob::Float64
  mild_detection_delay::Float64

  backward_tracing_prob::Float32
  backward_detection_delay::TimeDiff

  forward_tracing_prob::Float32
  forward_detection_delay::TimeDiff

  quarantine_length::Float32
  testing_time::TimeDiff

  phone_tracing_params::Union{Nothing, PhoneTracingParams}

  infection_modulation::Union{Nothing, InfectionModulation}
  mild_detection_modulation::Union{Nothing, InfectionModulation}
  forward_tracing_modulation::Union{Nothing, InfectionModulation}
  backward_tracing_modulation::Union{Nothing, InfectionModulation}

  screening_params::Union{Nothing, ScreeningParams}
  spreading_params::Union{Nothing, SpreadingParams}

  household_params::Union{Nothing, HouseholdParams}
end

numindividuals(params::SimParams) = length(params.household_ptrs)

straininfectivity(params::SimParams, strain::StrainKind) = straininfectivity(params.strain_infectivity_table, strain)
isimmune(state::SimState, subject_id::Integer, immunity::ImmunityState) = immunited(immunity) && time(state) <= state.individuals[subject_id].infections_immuity

householdof(params::SimParams, person_id::Integer) = UnitRange(params.household_ptrs[person_id]...)
age(params::SimParams, person_id::Integer) = params.ages[person_id]
gender(params::SimParams, person_id::Integer) = params.genders[person_id]

socialcompetence(params::SimParams, person_id::Integer) =
  nothing!=params.friendship_kernel_params && socialcompetence(params.friendship_kernel_params, person_id)

ishealthcare(params::SimParams, person_id::Integer) =
  nothing!=params.hospital_kernel_params && ishealthcare(params.hospital_kernel_params, person_id)
uses_phone_tracing(params::SimParams, person_id::Integer) =
  nothing!=params.phone_tracing_params && uses_phone_tracing(params.phone_tracing_params, person_id)

spreading(params::SimParams, person_id::Integer) = isnothing(params.spreading_params) ? 1.0 : spreading(params.spreading_params, person_id)

milddetectionprob(s::SimState, p::SimParams) = p.mild_detection_prob * evalmodulation(p.mild_detection_modulation, s, p)
forwardtracingprob(s::SimState, p::SimParams) = p.forward_tracing_prob * evalmodulation(p.forward_tracing_modulation, s, p)
backwardtracingprob(s::SimState, p::SimParams) = p.backward_tracing_prob * evalmodulation(p.backward_tracing_modulation, s, p)

milddetectiondelaydist(p::SimParams) = Exponential(p.mild_detection_delay)
forwarddetectiondelaydist(params::SimParams) = Exponential(params.forward_detection_delay)
backwarddetectiondelaydist(params::SimParams) = Exponential(params.backward_detection_delay)

function load_params(
  rng=MersenneTwister(0);
  population::DataFrame,

  infection_modulation_name::Union{Nothing,AbstractString} = nothing,
  infection_modulation_params::NamedTuple = NamedTuple{}(),

  mild_detection_modulation_name::Union{Nothing, AbstractString} = nothing,
  mild_detection_modulation_params::NamedTuple = NamedTuple{}(),

  forward_tracing_modulation_name::Union{Nothing, AbstractString} = nothing,
  forward_tracing_modulation_params::NamedTuple = NamedTuple{}(),

  backward_tracing_modulation_name::Union{Nothing, AbstractString} = nothing,
  backward_tracing_modulation_params::NamedTuple = NamedTuple{}(),

  hospitalization_time_ratio::Float64=1.0,
  hospitalization_multiplier::Float64=1.0, 
  death_multiplier::Float64=1.0,

  kwargs...
  )

  individuals_df::DataFrame = population

  progression_params::ProgressionParams = make_progression_params(hospitalization_time_ratio, hospitalization_multiplier, death_multiplier)

  infection_modulation = make_infection_modulation(infection_modulation_name; infection_modulation_params...)
  mild_detection_modulation = make_infection_modulation(mild_detection_modulation_name; mild_detection_modulation_params...)
  forward_tracing_modulation = make_infection_modulation(forward_tracing_modulation_name; forward_tracing_modulation_params...)
  backward_tracing_modulation = make_infection_modulation(backward_tracing_modulation_name; backward_tracing_modulation_params...)

  make_params(
    rng;
    individuals_df = individuals_df,
    progression_params = progression_params,
    infection_modulation = infection_modulation,
    mild_detection_modulation = mild_detection_modulation,
    forward_tracing_modulation = forward_tracing_modulation,
    backward_tracing_modulation = backward_tracing_modulation,
    kwargs...
  )
end

function make_params(
  rng::AbstractRNG=MersenneTwister(0);
  individuals_df::DataFrame,
  progression_params::ProgressionParams,

  infection_modulation=nothing,
  mild_detection_modulation=nothing,
  forward_tracing_modulation=nothing,
  backward_tracing_modulation=nothing,

  constant_kernel_param::Float64=1.0,
  household_kernel_param::Float64=1.0,

  hospital_detections::Bool=true,
  mild_detection_prob::Float64=0.0,
  mild_detection_delay::Float64=2.0,

  backward_tracing_prob::Float64=0.0,
  backward_detection_delay::Float64=1.0,

  forward_tracing_prob::Float64=0.0,
  forward_detection_delay::Float64=1.0,

  testing_time::Float64=1.0,
  quarantine_length::Float64=14.0,

  age_coupling_param::Union{Nothing, Real}=nothing,
  age_coupling_thresholds::Union{Nothing, AbstractArray{T} where T<:Real}=nothing,
  age_coupling_weights::Union{Nothing, AbstractMatrix{T} where T<:Real}=nothing,
  age_coupling_use_genders::Bool=false,

  screening_params::Union{Nothing,ScreeningParams}=nothing,
  household_params::Union{Nothing,HouseholdParams}=nothing,

  spreading_alpha::Union{Nothing,Real}=nothing,
  spreading_x0::Real=1,
  spreading_truncation::Real=Inf,

  phone_tracing_usage::Real=0.0,
  phone_detection_delay::Real=0.25,
  phone_tracing_usage_by_household::Bool=false,

  british_strain_multiplier::Real=1.70,
  delta_strain_multiplier::Real=1.7*1.5,
  omicron_strain_multiplier::Real=1.7*1.5*2.0,

  hospital_kernel_param::Float64=0.0,
  healthcare_detection_prob::Float64=0.8,
  healthcare_detection_delay::Float64=1.0,
)
  sort!(individuals_df, :household_index)

  num_individuals = individuals_df |> nrow

  household_ptrs = make_household_ptrs(individuals_df.household_index)

  strain_infectivity_table = make_infectivity_table(british_multiplier=british_strain_multiplier, delta_multiplier=delta_strain_multiplier,omicron_multiplier=omicron_strain_multiplier)

  age_coupling_kernel_params =
    if nothing === age_coupling_weights && nothing === age_coupling_thresholds && nothing === age_coupling_param; nothing
    elseif age_coupling_weights !== nothing && age_coupling_thresholds !== nothing
      @assert minimum(individuals_df.age) >= 0
      @assert maximum(individuals_df.age) * (age_coupling_use_genders+1) < typemax(GroupIdx)
      AgeCouplingParams(
        individuals_df.age,
        age_coupling_use_genders === nothing ? individuals_df.gender : nothing,
        age_coupling_thresholds, age_coupling_weights,
        age_coupling_param
        )
    else error("age couplig params not fully given")
    end

  hospital_kernel_params =
    if 0 == hospital_kernel_param; nothing
    elseif 0.0 < hospital_kernel_param;
      HospitalInfectionParams(
        individuals_df.ishealthcare,
        (UInt32(1):UInt32(num_individuals))[individuals_df.ishealthcare],
        hospital_kernel_param,
        healthcare_detection_prob,
        healthcare_detection_delay
      )
    else error("hospital_kernel_param must be postive or 0, got $hospital_kernel_param")
    end

  phone_tracing_params =
    if 0 == phone_tracing_usage; nothing
    elseif 0.0 < phone_tracing_usage <= 1.0
      phone_tracing_usage_by_household ?
        PhoneTracingParams(rng, num_individuals, phone_tracing_usage, phone_detection_delay, 1, household_ptrs) :
        PhoneTracingParams(rng, num_individuals, phone_tracing_usage, phone_detection_delay, 1)
    else error("tracing_app_usage must be nonnegative, got $phone_tracing_usage")
    end

  spreading_params =
    if nothing === spreading_alpha; nothing
    elseif 0.0 < spreading_alpha
      SpreadingParams(rng, num_individuals, alpha=spreading_alpha, x0=spreading_x0, truncation=spreading_truncation)
    else error("spreading_alpha must be larger than 0, got $spreading_alpha")
    end

  household_params =
    if nothing === household_params; HouseholdParams()
    elseif 0.0 <= household_params.quarantine_prob && 1.0 >= household_params.quarantine_prob && 0.0 <= household_params.trace_prob && 1.0 >= household_params.trace_prob
      household_params
    else error("household_params has invalid values for probas")
    end

  params = SimParams(
    household_ptrs,
    individuals_df.age,
    individuals_df.gender,
    reinterpret(ImmunityRand, rand(rng, UInt32, num_individuals)),

    progression_params,
    strain_infectivity_table,

    constant_kernel_param,
    household_kernel_param,
    age_coupling_kernel_params,
    hospital_kernel_params,

    hospital_detections,
    mild_detection_prob,
    mild_detection_delay,

    backward_tracing_prob,
    backward_detection_delay,

    forward_tracing_prob,
    forward_detection_delay,

    quarantine_length, # quarantine length
    testing_time, # testing time

    phone_tracing_params,

    infection_modulation,
    mild_detection_modulation,
    forward_tracing_modulation,
    backward_tracing_modulation,

    screening_params,
    spreading_params,
    household_params
  )
  params
end

function saveparams(dict, p::SimParams)
  dict["num_individuals"] = numindividuals(p)
  dict["constant/kernel_param"] = p.constant_kernel_param
  dict["household/kernel_param"] = p.household_kernel_param
  dict["quarantine/duration"] = p.quarantine_length
  dict["detections/hospital"] = p.hospital_detections
  dict["detections/mild"] = p.mild_detection_prob

  dict["infection_modulation"] = p.infection_modulation

  dict["tracing/backward_prob"] = p.backward_tracing_prob
  dict["tracing/backward_detection_delay"] = p.backward_detection_delay
  dict["tracing/forward_tracing_prob"] = p.forward_tracing_prob
  dict["tracing/forward_detection_delay"] = p.forward_detection_delay
  dict["tracing/testing_time"] = p.testing_time

  saveparams(dict, p.progressions, "progressions/")
  nothing===p.hospital_kernel_params || saveparams(dict, p.hospital_kernel_params, "hospital/")
  nothing===p.phone_tracing_params || saveparams(dict, p.phone_tracing_params, "phone_tracing/")
  nothing===p.spreading_params || saveparams(dict, p.spreading_params, "spreading")
  nothing
end
