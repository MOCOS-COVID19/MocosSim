struct Progression
  severity::Severity
  # times given with respect to the infection time
  incubation_time::Union{Missing, TimeDiff}
  mild_symptoms_time::Union{Missing,TimeDiff}
  severe_symptoms_time::Union{Missing,TimeDiff}
  #critical_symptoms_time::Float32
  recovery_time::Union{Missing,TimeDiff}
  death_time::Union{Missing, TimeDiff}
  #Progression(severity::Severity, incubation_time::Real, mild_time::Real, severe_time::Real, recovery_time) = incubation < mild_time < severe_time < recovery_time
end

Progression() = Progression(UndefinedSeverity, missing, missing, missing, missing, missing)