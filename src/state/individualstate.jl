struct IndividualState
    health::HealthState
    freedom::FreedomState
    detected::DetectionStatus
    immunity::ImmunityState
    quarantine_level::SafeUInt16
  end

  IndividualState() = IndividualState(
    Healthy,
    Free,
    Undetected,
    NoImmunity,
    0
  )

  show(io::IO, s::IndividualState) = print(io, "(",s.health, ", ", s.freedom, ", ", s.detected, ", ", s.immunity, ", ", s.quarantine_level, ")")
