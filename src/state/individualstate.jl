struct IndividualState
    health::HealthState
    freedom::FreedomState
    detected::DetectionStatus
    immunity::ImmunityState
    immunization_day::TimeDay
    quarantine_level::SafeUInt16
    strain::StrainKind
    infections_immuity::TimePoint
    severe_immunity::TimePoint
  end

  IndividualState() = IndividualState(
    Healthy,
    Free,
    Undetected,
    NoImmunity,
    0,
    0,
    NullStrain,
    0,
    0
  )

  show(io::IO, s::IndividualState) = print(io, "(",s.health, ", ", s.freedom, ", ", s.detected, ", ", s.immunity, ", ", s.imunization_day, ", ", s.quarantine_level,  ", ", s.strain,  ")")
