@enum EventKind::UInt8 begin

# the underlying values are also priorites in case the times are equal
# therefore the order of definition implies priority

  QuarantinedEvent
  DetectionEvent
  #DetectedFromQuarantineEvent
  HomeTreatmentSuccessEvent
  QuarantineEndEvent

  GoHospitalEvent

  #DetectedOutsideQuarantineEvent
  #DetectedFromTracingEvent

  TracedEvent
  ReleasedEvent

  # the progression events have low priority to let the immediate actions execute
  BecomeInfectiousEvent

  OutsideInfectionEvent
  TransmissionEvent
  OutsideTransmissionEvent

  MildSymptomsEvent
  HomeTreatmentEvent

  SevereSymptomsEvent
  CriticalSymptomsEvent
  RecoveryEvent
  DeathEvent

  ScreenigEvent

  InvalidEvent # should not be executed

end

struct Event
  time::TimePoint                 # 4 bytes
  subject_id::PersonIdx           # 4 bytes
  source_id::PersonIdx            # 4 bytes
  event_kind::EventKind           # 1 byte
  extra::UInt8
  strain::StrainKind

  #contact_kind::ContactKind       # 1 byte
  #extension::Bool                 # 1 byte
  #detection_kind::DetectionKind   # 1 byte
                                  # alignment = 4 bytes

  Event() = new(0.0, 0, 0, InvalidEvent, 0, NullStrain)
  Event(::Val{E}, time::Real, subject::Integer) where E = new(time, subject, 0, E, 0, NullStrain)
  Event(::Val{OutsideInfectionEvent}, time::Real, subject::Integer, strain::StrainKind) = new(time, subject, 0, OutsideInfectionEvent, UInt8(OutsideContact), strain)
  Event(::Val{OutsideTransmissionEvent}, time::Real, subject::Integer, strain::StrainKind) = new(time, subject, 0, OutsideTransmissionEvent, UInt8(OutsideContact), strain)
  Event(::Val{TransmissionEvent}, ::Real, ::Integer) = error("source and contact kind needed for transmission event")
  Event(::Val{TransmissionEvent}, time::Real, subject::Integer, source::Integer, contact_kind::ContactKind, strain::StrainKind) = new(time, subject, source, TransmissionEvent, UInt8(contact_kind), strain)
  Event(::Val{QuarantinedEvent}, time::Real, subject::Integer, extension::Bool) = new(time, subject, 0, QuarantinedEvent, UInt8(extension), NullStrain)
  Event(::Val{TracedEvent}, ::Real, ::Integer) = error("source and tracing kind must be given for TracedEvent")
  Event(::Val{TracedEvent}, time::Real, subject::Integer, source::Integer, tracing_kind::TracingKind) = new(time, subject, source, TracedEvent, UInt8(tracing_kind), NullStrain)
  Event(::Val{DetectionEvent}, ::Real, ::Integer) = error("detection kind must be given for detection event")
  Event(::Val{DetectionEvent}, time::Real, subject::Integer, detectionkind::DetectionKind) = new(time, subject, 0, DetectionEvent, UInt8(detectionkind), NullStrain)
  Event(::Val{ScreenigEvent}, time::Real) = new(time, 0, 0, ScreenigEvent, 0, NullStrain)
end

time(event::Event) = event.time
subject(event::Event) = event.subject_id
source(event::Event) = event.source_id
kind(event::Event) = event.event_kind

contactkind(event::Event) = istransmission(event) ? ContactKind(event.extra) : NoContact
detectionkind(event::Event) = isdetection(event) ? DetectionKind(event.extra) : NoDetection
extension(event::Event) = isquarantine(event) ? Bool(event.extra) : false
tracingkind(event::Event) = istracing(event) ? TracingKind(event.extra) : NotTraced

strainkind(event::Event) = istransmission(event) ? event.strain : NullStrain

function show(io::IO, event::Event)
  print(io, time(event), ":", kind(event), " ", subject(event))
  if TransmissionEvent == kind(event) || OutsideInfectionEvent == kind(event)
    print(io, " <= ", source(event), " ", contactkind(event), " ", strainkind(event))
  elseif QuarantinedEvent == kind(event)
    print(io, " extension=", extension(event))
  elseif TracedEvent == kind(event)
    print(io, " <= ", source(event))
  elseif InvalidEvent == kind(event)
    print(io, " <= ", source(event), " ", kind(event), " ", event.extra)
  end
end

isdetection(ek::EventKind) = ek == DetectionEvent
isdetection(e::Event) = isdetection(kind(e))
istransmission(ek::EventKind) = ek == TransmissionEvent || ek == OutsideInfectionEvent
istransmission(e::Event) = istransmission(kind(e))
isquarantine(ek::EventKind) = ek == QuarantinedEvent || ek == QuarantineEndEvent
isquarantine(e::Event) = isquarantine(kind(e))
istracing(ek::EventKind) = ek == TracedEvent
istracing(e::Event) = istracing(kind(e))
ishospitalization(ek::EventKind) = ek == GoHospitalEvent
ishospitalization(e::Event) = ishospitalization(kind(e))
isdeath(ek::EventKind) = ek == DeathEvent
isdeath(e::Event) = isdeath(kind(e))