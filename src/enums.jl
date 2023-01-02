@enum Severity::UInt8 begin
  UndefinedSeverity=0
  Asymptomatic=1
  Mild
  Severe
  Critical
end

@enum HealthState::UInt8 begin # 3 bits
  Healthy
  Incubating
  Infectious
  MildSymptoms
  SevereSymptoms
  CriticalSymptoms
  # Recovered
  Dead
end

@enum FreedomState::UInt8 begin # 3 bits
  Free
  HomeQuarantine
  HomeTreatment
  Hospitalized
  # Released
end

@enum DetectionStatus::UInt8 begin #2 bits
  Undetected
  UnderObservation
  TestPending
  Detected
end

@enum ContactKind::UInt8 begin  # 3 bits
  NoContact=0
  HouseholdContact
  HospitalContact
  AgeCouplingContact
  ConstantKernelContact
  OutsideContact
end

@enum DetectionKind::UInt8 begin
  NoDetection=0
  OutsideQuarantineDetection=1
  FromQuarantineDetection
  FromTracingDetection
end

@enum TracingKind::UInt8 begin
  NotTraced = 0
  MildCaseTraced #unused yet
  HospitalTraced #unused yet
  QuarantineTraced #unused yet
  ClassicalTraced
  PhoneTraced
end

@enum StrainKind::UInt8 begin
  NullStrain = 0
  ChineseStrain
  BritishStrain
  DeltaStrain
  OmicronStrain
end

const NUM_STRAINS = length(instances(StrainKind)) - 1

@enum ImmunityState::UInt8 begin
  NullImmunity = 0
  NoImmunity
  NaturalImmunity
  VecVacImmunity
  MRNAVacImmunity
end

const NUM_IMMUNITIES = length(instances(ImmunityState)) - 1
