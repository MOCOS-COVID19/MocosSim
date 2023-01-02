struct IndividualState
  health::HealthState
  freedom::FreedomState
  detected::DetectionStatus
  immunity::ImmunityState
  immunization_day::TimeDay
  quarantine_end::TimeDay     # the first day quarantine-free, 0 for no quarantine
  #quarantine_level::SafeUInt16
end

  IndividualState() = IndividualState(
    Healthy,
    Free,
    Undetected,
    NoImmunity,
    0,
    0
  )

  show(io::IO, s::IndividualState) = print(io, "(",s.health, ", ", s.freedom, ", ", s.detected, ", ", s.immunity, ", ", s.immunization_day, ", ", s.quarantine_end, ")")
