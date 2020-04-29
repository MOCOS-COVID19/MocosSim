struct Progression
    severity::Severity
    # times given with respect to the infection time
    incubation_time::TimeDiff
    mild_symptoms_time::Union{Missing,TimeDiff}
    severe_symptoms_time::Union{Missing,TimeDiff}
    #critical_symptoms_time::Float32
    recovery_time::Union{Missing,TimeDiff}
    death_time::Union{Missing, TimeDiff}
    #Progression(severity::Severity, incubation_time::Real, mild_time::Real, severe_time::Real, recovery_time) = incubation < mild_time < severe_time < recovery_time
end
Progression() = Progression(Asymptomatic, missing, missing, missing, missing, missing)

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

function sample_severity(rng::AbstractRNG, age::Real)
  dist = severity_dists[1]
  if age < 0;       error("age should be non-negative")
  elseif age < 40;  dist = severity_dists[1]
  elseif age < 50;  dist = severity_dists[2]    
  elseif age < 60;  dist = severity_dists[3]
  elseif age < 70;  dist = severity_dists[4]
  elseif age < 80;  dist = severity_dists[5]
  else;             dist = severity_dists[6]
  end
  severity_int = rand(rng, dist)
  severity = rand(rng, dist) |> Severity
end

function sample_if_death(rng::AbstractRNG, severity::Severity)
  severity_int = severity |> UInt8
  death_prob = death_probs[severity_int]
  rand(rng) < death_prob
end

@inline function sample_progression(rng::AbstractRNG, age::Real, 
    dist_incubation, 
    dist_symptom_onset, 
    dist_hospitalization,         
    dist_mild_recovery,
    dist_severe_recovery,
    dist_death_time
  )

  severity = sample_severity(rng, age)
 
  mild_symptoms_time = missing
  severe_symptoms_time = missing
  recovery_time = missing
  death_time = missing
    
  incubation_time = rand(rng, dist_incubation)
   
  mild_symptoms_time = incubation_time + rand(rng, dist_symptom_onset)
  
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
    recovery_time = mild_symptoms_time + rand(rng, dist_mild_recovery)
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

function resample_progressions!(rng::AbstractRNG, 
                                progressions, ages::AbstractArray{T} where T <: Real,       
                                dist_incubation_time, 
                                dist_symptom_onset_time, 
                                dist_hospitalization_time,
                                dist_mild_recovery_time,
                                dist_severe_recovery_time,
                                dist_death_time
                                )

  for i in 1:length(ages)
    progressions[i] = Simulation.sample_progression(rng, ages[i],         
      dist_incubation_time, 
      dist_symptom_onset_time, 
      dist_hospitalization_time,
      dist_mild_recovery_time,
      dist_severe_recovery_time,
      dist_death_time)
  end
  progressions
end