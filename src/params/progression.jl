struct ProgressionParams
  # the distribution are subject to change
  dist_incubation_time::LogNormal{Float64}
  dist_symptom_onset_time::Gamma{Float64}
  dist_hospitalization_time::Gamma{Float64}
  dist_mild_recovery_time::Uniform{Float64}
  dist_severe_recovery_time::Uniform{Float64}
  dist_death_time::LogNormal{Float64}
end

function make_progression_params()
  #this data needs updating
  dist_incubation_time = LogNormal(1.3669786931887833, 0.5045104580676582)
  dist_symptom_onset_time = Gamma(0.8738003969079596, 2.9148873266517685)
  dist_hospitalization_time = Gamma(1.1765988120148885, 2.6664347368236787)
  dist_mild_recovery_time = Uniform(11, 17)
  dist_severe_recovery_time = Uniform(4*7, 8*7)
  dist_death_time = LogNormal(2.610727719719777, 0.44476420066780653)

  ProgressionParams(
    dist_incubation_time,
    dist_symptom_onset_time,
    dist_hospitalization_time,
    dist_mild_recovery_time,
    dist_severe_recovery_time,
    dist_death_time
  )
end

# Asymptotic not allowed as yet (some asserts prevent them)
const severity_dists = [
  Categorical([0,  0.852,  0.144,  0.004]),
  Categorical([0,  0.848,  0.144,  0.008]),
  Categorical([0,  0.832,  0.141,  0.027]),
  Categorical([0,  0.793,  0.134,  0.073]),
  Categorical([0,  0.716,  0.121,  0.163]),
  Categorical([0,  0.596,  0.102,  0.302]),
]

const death_probs = [0.0, 0.0, 0.0, 0.49]

const vaccination_effectiveness = [0.0, 0.75, 0.875, 0.92]

function sample_severity(rng::AbstractRNG, age::Real, immunity::ImmunityState)
  dist = severity_dists[1]
  if age < 0;       error("age should be non-negative")
  elseif age < 40;  dist = severity_dists[1]
  elseif age < 50;  dist = severity_dists[2]
  elseif age < 60;  dist = severity_dists[3]
  elseif age < 70;  dist = severity_dists[4]
  elseif age < 80;  dist = severity_dists[5]
  else;             dist = severity_dists[6]
  end

  severity = rand(rng, dist) |> Severity
  severity_int = severity |> UInt8
  @assert length(vaccination_effectiveness) == length(instances(Severity)) - 1
  #reduction severity for immunited subject with some probability
  if severity > Asymptomatic && immunited(immunity) && rand(rng) < vaccination_effectiveness[severity_int]
    severity_int = severity_int - 0x01
    severity = severity_int |> Severity
  end
  severity
end

function sample_if_death(rng::AbstractRNG, severity::Severity)
  severity_int = severity |> UInt8
  death_prob = death_probs[severity_int]
  rand(rng) < death_prob
end

function sample_progression(rng::AbstractRNG, progression_data::ProgressionParams, age::Real, gender::Bool, immunity::ImmunityState, time_since_immunization::Real, strain::StrainKind)
  #TODO expecting updated progression generation soon
  sample_progression(
    rng,
    age,
    immunity,
    progression_data.dist_incubation_time,
    progression_data.dist_symptom_onset_time,
    progression_data.dist_hospitalization_time,
    progression_data.dist_mild_recovery_time,
    progression_data.dist_severe_recovery_time,
    progression_data.dist_death_time
  )
end


@inline function sample_progression(rng::AbstractRNG, age::Real, immunity::ImmunityState,
    dist_incubation,
    dist_symptom_onset,
    dist_hospitalization,
    dist_mild_recovery,
    dist_severe_recovery,
    dist_death_time
  )

  severity = sample_severity(rng, age, immunity)

  mild_symptoms_time = missing
  severe_symptoms_time = missing
  recovery_time = missing
  death_time = missing

  incubation_time = rand(rng, dist_incubation)

  if (severity==Mild) || (severity==Severe) || (severity==Critical)
    mild_symptoms_time = incubation_time + rand(rng, dist_symptom_onset)
  end

  severe_symptoms_time = missing
  recovery_time = missing
  death_time = missing
  if (severity==Severe) || (severity==Critical)
    severe_symptoms_time = incubation_time + rand(rng, dist_hospitalization)
    if severe_symptoms_time <= mild_symptoms_time
      mild_symptoms_time = missing
    end
    recovery_time = severe_symptoms_time + rand(rng, dist_severe_recovery)
  else
    if (severity==Mild)
      recovery_time = mild_symptoms_time + rand(rng, dist_mild_recovery)
    else
      if immunited(immunity) # && (severity==Asymptomatic)
        recovery_time = incubation_time  # TODO if we want this # rand(rng, dist_mild_recovery)
      else # now only asymptomatic, but not vaccinated
        recovery_time = rand(rng, dist_mild_recovery)
      end
    end
  end

  if sample_if_death(rng, severity)
    death_time = incubation_time + rand(rng, dist_death_time)
    if death_time < severe_symptoms_time
      death_time = severe_symptoms_time 
    end
    recovery_time = missing
  end

  Progression(
    severity,
    incubation_time |> TimeDiff,
    ismissing(mild_symptoms_time) ? missing : TimeDiff(mild_symptoms_time),
    ismissing(severe_symptoms_time) ? missing : TimeDiff(severe_symptoms_time),
    ismissing(recovery_time) ? missing : TimeDiff(recovery_time),
    ismissing(death_time) ? missing : TimeDiff(death_time)
  )
end

function resample!(
  rng::AbstractRNG,
  progressions, ages::AbstractArray{T} where T <: Real,
  dist_incubation_time,
  dist_symptom_onset_time,
  dist_hospitalization_time,
  dist_mild_recovery_time,
  dist_severe_recovery_time,
  dist_death_time)

  for i in 1:length(ages)
    progressions[i] = sample_progression(rng, ages[i],
      dist_incubation_time,
      dist_symptom_onset_time,
      dist_hospitalization_time,
      dist_mild_recovery_time,
      dist_severe_recovery_time,
      dist_death_time)
  end
  progressions
end

storagefloat(t::Union{Missing,TimeDiff})::Float32 = ismissing(t) ? NaN32 : Float32(t)
storagefloat(p::Progression, symbol::Symbol) = getproperty(p, symbol) |> storagefloat

function saveparams(dict, progressions::AbstractVector{Progression}, prefix::AbstractString="")
    dict[prefix*"severities"] = getproperty.(progressions, :severity)
    dict[prefix*"incubation_times"] = storagefloat.(progressions, :incubation_time)
    dict[prefix*"mild_symptoms_times"] = storagefloat.(progressions, :mild_symptoms_time)
    dict[prefix*"severe_symptoms_times"] = storagefloat.(progressions, :severe_symptoms_time)
    dict[prefix*"recovery_times"] = storagefloat.(progressions, :recovery_time)
    dict[prefix*"death_times"] = storagefloat.(progressions, :death_time)
end