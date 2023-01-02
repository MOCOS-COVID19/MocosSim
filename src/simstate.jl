include("state/individualstate.jl")
include("state/progression.jl")
include("state/runningstats.jl")

abstract type AbstractSimState end

# the mutable part of the simulation
mutable struct SimState <: AbstractSimState
  rng::MersenneTwister
  time::TimePoint
  queue::EventQueue
  individuals::Vector{IndividualState}
  progressions::Vector{Progression}
  forest::InfectionForest
  stats::RunningStats

  SimState(rng::AbstractRNG, num_individuals::Integer) = num_individuals <= 0 || num_individuals > typemax(UInt32) ? error("number of individuals must be positive and smaller than $(typemax(UInt32))") : 
    new(
      rng,

      0.0,
      EventQueue(),

      fill(IndividualState(), num_individuals),
      fill(Progression(), num_individuals),

      InfectionForest(num_individuals),

      RunningStats()
    )
end

SimState(num_individuals::Integer; seed::Integer=0) = SimState(MersenneTwister(seed), num_individuals)

function show(io::IO, state::SimState)
  print(io, "Simulation state for ", length(state.individuals), " individuals")
end

function reset!(state::SimState, rng::AbstractRNG)
  state.rng = rng
  state.time=0
  empty!(state.queue)
  reset!(state.forest)
  fill!(state.individuals, IndividualState())
  fill!(state.progressions, Progression())
  reset!(state.stats)
  state
end

time(state::SimState) = state.time
numindividuals(state::SimState) = length(state.individuals)

reset!(state::SimState) = reset!(state::SimState, state.rng)
reset!(state::SimState, seed::Integer) = reset!(state, MersenneTwister(seed))

individualstate(state::SimState, person_id::Integer) = state.individuals[person_id]
health(state::SimState, person_id::Integer)::HealthState = state.individuals[person_id].health
freedom(state::SimState, person_id::Integer)::FreedomState = state.individuals[person_id].freedom

detected(state::SimState, person_id::Integer)::DetectionStatus = state.individuals[person_id].detected
isdetected(state::SimState, person_id::Integer)::Bool = (Detected == detected(state, person_id))

subjecthealth(state::SimState, event::Event)::HealthState = health(state, subject(event))
subjectfreedom(state::SimState, event::Event)::FreedomState = freedom(state, subject(event))

sourcehealth(state::SimState, event::Event)::HealthState = health(state, source(event))
sourcefreedom(state::SimState, event::Event)::FreedomState = freedom(state, source(event))

recentstrainof(state::SimState, person_id::Integer) = recentstrainof(state.forest, person_id)
immunityof(state::SimState, person_id::Integer)::ImmunityState = state.individuals[person_id].immunity
immunizationday(state::SimState, person_id::Integer) = state.individuals[person_id].immunization_day
timesinceimmunization(state::SimState, person_id::Integer)::TimePoint = time(state) - TimePoint(immunizationday(state, person_id))

function progressionof(state::SimState, person_id::Integer)::Progression
  progression = state.progressions[person_id]
  @assert progression.severity !== UndefinedSeverity
  progression
end
numdetected(state::SimState) = numdetected(state.stats)
numdead(state::SimState) = numdead(state.stats)

recentforwardinfectionsof(state::SimState, person_id::Integer) = recentforwardinfectionsof(state.forest, person_id)
recentbackwardinfection(state::SimState, person_id::Integer) = recentbackwardinfectionof(state.forest, person_id)

function sethealth!(state::SimState, person_id::Integer, new_health::HealthState)
  orig = state.individuals[person_id]
  @assert orig.health <= new_health || Healthy == new_health
  state.individuals[person_id] = @set orig.health = new_health
  nothing
end

function setfreedom!(state::SimState, person_id::Integer, new_freedom::FreedomState)
  orig = state.individuals[person_id]
  state.individuals[person_id] = @set orig.freedom = new_freedom
  nothing
end

function setdetected!(state::SimState, person_id::Integer, new_detected::DetectionStatus)
  orig = state.individuals[person_id]
  @assert orig.detected <= new_detected
  state.individuals[person_id] = @set orig.detected = new_detected
  nothing
end

function setimmunity!(state::SimState, person_id::Integer, new_immunity::ImmunityState, time::Real)
  orig = state.individuals[person_id]
  state.individuals[person_id] = @set orig.immunity = new_immunity
end

setimmunity!(state::SimState, person_id::Integer, new_immunity::ImmunityState) =
  setimmunity!(state, person_id, new_immunity, time(state))

function clearprogression!(state::SimState, person_id::Integer)
  @assert state.progressions[person_id].severity != UndefinedSeverity
  state.progressions[person_id] = Progression()
  nothing
end

function setprogression!(state::SimState, person_id::Integer, progression::Progression)
  @assert state.progressions[person_id].severity == UndefinedSeverity
  state.progressions[person_id] = progression
  nothing
end

#quarantine_level(state::SimState, person_id::Integer) = state.individuals[person_id].quarantine_level
#isquarantined(state::SimState, person_id::Integer)::Bool = quarantine_level(state, person_id) != 0

#function quarantine_advance!(state::SimState, person_id::Integer, adv_val::Integer)
#  orig = state.individuals[person_id]
#  if adv_val >= 0
#    state.individuals[person_id] = @set orig.quarantine_level = orig.quarantine_level+adv_val
#  else
#    state.individuals[person_id] = @set orig.quarantine_level = orig.quarantine_level-(-adv_val)
#  end
#  nothing
#end
quarantine_end(state::SimState, person_id::Integer) = TimePoint(state.individuals[person_id].quarantine_end)
isquarantined(state::SimState, person_id::Integer)::Bool = time(state) < quarantine_end(state, person_id)

function quarantine_prolong!(state::SimState, person_id::Integer, new_time::Real)
  orig = state.individuals[person_id]
  new_time = ceil(TimeDay, new_time)
  @assert orig.quarantine_end <= new_time
  state.individuals[person_id] = @set orig.quarantine_end = new_time;
  nothing
end

function quarantine_cancel!(state::SimState, person_id::Integer)
  orig = state.individuals[person_id]
  #state.individuals[person_id] = @set orig.quarantine_level = 0
  state.individuals[person_id] = @set orig.quarantine_end = 0
  nothing
end

setsubjecthealth!(state::SimState, event::Event, health::HealthState) = sethealth!(state, subject(event), health)
setsubjectfreedom!(state::SimState, event::Event, freedom::FreedomState) = setfreedom!(state, subject(event), freedom)

registerinfection!(state::SimState, infection::Event) = push!(state.forest, infection)

function initialfeed!(state::SimState, num_initial_infections::Integer, strain::StrainKind=ChineseStrain)
  @assert 0 == time(state)

  N = length(state.individuals)
  individuals = 1:N

  for _ in 1:num_initial_infections
    person_id = sample(state.rng, individuals)
    event = Event(Val(OutsideInfectionEvent), 0.0, person_id, strain)
    push!(state.queue, event)
  end

end

saveparams(dict, state::SimState, prefix::AbstractString="") = saveparams(dict, state.forest, prefix)