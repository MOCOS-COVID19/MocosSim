struct IndividualState
    health::HealthState
    freedom::FreedomState
    detected::DetectionStatus
    infection_immunity::Bool
    severe_immunity::Bool
    immunization_day::TimeDay
    quarantine_level::SafeUInt16
    strain::StrainKind
  end

  IndividualState() = IndividualState(
    Healthy,
    Free,
    Undetected,
    false,
    false,
    0,
    0,
    NullStrain
  )

  show(io::IO, s::IndividualState) = print(io, "(",s.health, ", ", s.freedom, ", ", s.detected, ", ", s.immunity, ", ", s.imunization_day, ", ", s.quarantine_level,  ", ", s.strain,  ")")
