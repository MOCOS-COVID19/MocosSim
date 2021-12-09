@enum Severity::UInt8 UndefinedSeverity=0 Asymptomatic=1 Mild Severe Critical

@enum HealthState::UInt8 Healthy Incubating Infectious MildSymptoms SevereSymptoms CriticalSymptoms Recovered Dead # 3 bits

@enum FreedomState::UInt8 Free HomeQuarantine HomeTreatment Hospitalized Released # 3 bits

@enum DetectionStatus::UInt8 Undetected UnderObservation TestPending Detected #2 bits

@enum ContactKind::UInt8 NoContact=0 HouseholdContact HospitalContact AgeCouplingContact ConstantKernelContact OutsideContact # 3 bits

@enum DetectionKind::UInt8 NoDetection=0 OutsideQuarantineDetection=1 FromQuarantineDetection FromTracingDetection

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
