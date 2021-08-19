using Random
using Distributions
#using FunctionWrappers

const Age=UInt8

include("params/households.jl")
include("params/friendship.jl")
include("params/progression.jl")
include("params/hospital.jl")
include("params/phonetracing.jl")
include("params/spreading.jl")
include("params/strains.jl")

abstract type InfectionModulation end
abstract type ScreeningParam end

abstract type InfectionTravels end

struct SimParams
  household_ptrs::Vector{Tuple{PersonIdx,PersonIdx}}  # (i1,i2) where i1 and i2 are the indices of first and last member of the household

  ages::Vector{Age}
  genders::BitVector

  progressions::Vector{Progression} # not sure if progressions should be there

  strain_table::StrainTable

  hospital_kernel_params::Union{Nothing, HospitalInfectionParams}  # nothing if hospital kernel not active
  friendship_kernel_params::Union{Nothing, FriendshipKernelParams}

  hospital_detections::Bool
  mild_detection_prob::Float64

  backward_tracing_prob::Float32
  backward_detection_delay::TimeDiff

  forward_tracing_prob::Float32
  forward_detection_delay::TimeDiff

  quarantine_length::Float32
  testing_time::TimeDiff

  phone_tracing_params::Union{Nothing, PhoneTracingParams}

  infection_modulation_function::Union{Nothing,InfectionModulation}
  screening_params::Union{Nothing, ScreeningParam}
  infection_travels_function::Union{Nothing,InfectionTravels}
<<<<<<< HEAD
=======
  travels_frequency::TimePoint

>>>>>>> create event to outside trasmission
  spreading_params::Union{Nothing, SpreadingParams}
end

include("params/modulations.jl")
include("params/screening.jl")
include("params/travels.jl")

numindividuals(params::SimParams) = length(params.household_ptrs)
straindata(params::SimParams, strain::StrainKind) = getdata(params.strain_table, strain)

progressionof(params::SimParams, person_id::Integer) = params.progressions[person_id]
severityof(params::SimParams, person_id::Integer) = progressionof(params, person_id).severity
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

function load_params(rng=MersenneTwister(0);
        population::DataFrame,
        infection_modulation_name::Union{Nothing,AbstractString}=nothing,
        infection_modulation_params::NamedTuple=NamedTuple{}(),
        screening_params::Union{Nothing,ScreeningParams}=nothing,
        infection_travels_name::Union{Nothing,AbstractString}=nothing,
        infection_travels_params::NamedTuple=NamedTuple{}(),
        travels_frequency::TimePoint = 0.0,
        kwargs...
        )

  individuals_df::DataFrame = population

  num_individuals = individuals_df |> nrow

  dist_incubation_time = LogNormal(1.3669786931887833, 0.5045104580676582)
  dist_symptom_onset_time = Gamma(0.8738003969079596, 2.9148873266517685)
  dist_hospitalization_time = Gamma(1.1765988120148885, 2.6664347368236787)
  dist_mild_recovery_time = Uniform(11, 17)
  dist_severe_recovery_time = Uniform(4*7, 8*7)
  dist_death_time = LogNormal(2.610727719719777, 0.44476420066780653)

  progressions = Vector{Progression}(undef, num_individuals);
  resample!(rng, progressions, individuals_df.age,
    dist_incubation_time,
    dist_symptom_onset_time,
    dist_hospitalization_time,
    dist_mild_recovery_time,
    dist_severe_recovery_time,
    dist_death_time
  )

  infection_modulation_function = isnothing(infection_modulation_name) ? nothing : make_infection_modulation(infection_modulation_name; infection_modulation_params...)
  infection_travels_function = isnothing(infection_travels_name) ? nothing : make_infection_travels(infection_travels_name; infection_travels_params...)
  
  make_params(
    rng,
    individuals_df=individuals_df,
    progressions=progressions,
    infection_modulation_function=infection_modulation_function,
<<<<<<< HEAD
    screening_params=screening_params;
    infection_travels_function=infection_travels_function;
=======
    infection_travels_function=infection_travels_function,
    travels_frequency = travels_frequency;
>>>>>>> create event to outside trasmission
    kwargs...
  )
end

function make_params(
  rng::AbstractRNG=MersenneTwister(0);
  individuals_df::DataFrame,
  progressions::AbstractArray{Progression},

  infection_modulation_function=nothing,
  infection_travels_function=nothing,
  travels_frequency::TimePoint=0.0,

  screening_params::Union{Nothing,ScreeningParams}=nothing,

  hospital_kernel_param::Float64=0.0,
  healthcare_detection_prob::Float64=0.8,
  healthcare_detection_delay::Float64=1.0,

  constant_kernel_param::Float64=1.0,
  household_kernel_param::Float64=1.0,
  friendship_kernel_param::Float64=0.0,

  hospital_detections::Bool=true,
  mild_detection_prob::Float64=0.0,

  backward_tracing_prob::Float64=0.0,
  backward_detection_delay::Float64=1.0,

  forward_tracing_prob::Float64=0.0,
  forward_detection_delay::Float64=1.0,

  quarantine_length::Float64=14.0,
  testing_time::Float64=1.0,

  phone_tracing_usage::Real=0.0,
  phone_detection_delay::Real=0.25,
  phone_tracing_usage_by_household::Bool=false,

  spreading_alpha::Union{Nothing,Real}=nothing,
  spreading_x0::Real=1,
  spreading_truncation::Real=Inf,

  british_strain_multiplier::Real=1.70
)
  sort!(individuals_df, :household_index)

  num_individuals = individuals_df |> nrow

  @assert num_individuals == length(progressions)

  household_ptrs = make_household_ptrs(individuals_df.household_index)

  strain_table = make_strains(constant_kernel_param, household_kernel_param, british_multiplier=british_strain_multiplier)

  friendship_kernel_params =  if 0 == friendship_kernel_param; nothing
                        elseif 0.0 < friendship_kernel_param
                          FriendshipKernelParams(
                            friendship_kernel_param,
                            Age.(individuals_df.age),
                            individuals_df.gender,
                            individuals_df.social_competence)
                        else error("bad condition for friendship kernel")
                        end

  hospital_kernel_params =  if 0 == hospital_kernel_param; nothing
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

  phone_tracing_params = if 0 == phone_tracing_usage; nothing
                      elseif 0.0 < phone_tracing_usage <= 1.0
                        phone_tracing_usage_by_household ?
                          PhoneTracingParams(rng, num_individuals, phone_tracing_usage, phone_detection_delay, 1, household_ptrs) :
                          PhoneTracingParams(rng, num_individuals, phone_tracing_usage, phone_detection_delay, 1)
                      else error("tracing_app_usage must be nonnegative, got $phone_tracing_usage")
                      end

  spreading_params = if nothing == spreading_alpha; nothing
                  elseif 0.0 < spreading_alpha
                    SpreadingParams(rng, num_individuals, alpha=spreading_alpha, x0=spreading_x0, truncation=spreading_truncation)
                  else error("spreading_alpha must be larger than 0, got $spreading_alpha")
                  end

  params = SimParams(
    household_ptrs,
    individuals_df.age,
    individuals_df.gender,

    progressions,

    strain_table,

    friendship_kernel_params,
    hospital_kernel_params,

    hospital_detections,
    mild_detection_prob,

    backward_tracing_prob,
    backward_detection_delay,

    forward_tracing_prob,
    forward_detection_delay,

    quarantine_length, # quarantine length
    testing_time, # testing time

    phone_tracing_params,
    infection_modulation_function,
    screening_params,
    infection_travels_function,
    travels_frequency,
    spreading_params
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

  dict["infection_modulation"] = p.infection_modulation_function
  dict["infection_travels"] = p.infection_travels_function

  dict["tracing/backward_prob"] = p.backward_tracing_prob
  dict["tracing/backward_detection_delay"] = p.backward_detection_delay
  dict["tracing/forward_tracing_prob"] = p.forward_tracing_prob
  dict["tracing/forward_detection_delay"] = p.forward_detection_delay
  dict["tracing/testing_time"] = p.testing_time

  saveparams(dict, p.progressions, "progressions/")
  nothing==p.hospital_kernel_params || saveparams(dict, p.hospital_kernel_params, "hospital/")
  nothing==p.phone_tracing_params || saveparams(dict, p.phone_tracing_params, "phone_tracing/")
  nothing==p.spreading_params || saveparams(dict, p.spreading_params, "spreading")
  nothing
end
