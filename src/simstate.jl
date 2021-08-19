include("state/individualstate.jl")
include("state/runningstats.jl")

# the mutable part of the simulation
mutable struct SimState
  rng::MersenneTwister
  time::TimePoint
  queue::EventQueue
  individuals::Vector{IndividualState}
  forest::InfectionForest
  stats::RunningStats

  SimState(rng::AbstractRNG, num_individuals::Integer) = num_individuals <= 0 || num_individuals > typemax(UInt32) ? error("number of individuals must be positive and smaller than $(typemax(UInt32))") : 
    new(
      rng,

      0.0,
      EventQueue(),

      fill(IndividualState(), num_individuals),

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
  reset!(state.stats)
  state
end

time(state) = state.time
numindividuals(state::SimState) = length(state.individuals)

reset!(state::SimState) = reset!(state::SimState, state.rng)
reset!(state::SimState, seed::Integer) = reset!(state, MersenneTwister(seed))

individualstate(state::SimState, person_id::Integer) = state.individuals[person_id]
health(state::SimState, person_id::Integer)::HealthState = state.individuals[person_id].health
freedom(state::SimState, person_id::Integer)::FreedomState = state.individuals[person_id].freedom

quarantine_level(state::SimState, person_id::Integer) = state.individuals[person_id].quarantine_level
isquarantined(state::SimState, person_id::Integer)::Bool = quarantine_level(state, person_id) != 0
detected(state::SimState, person_id::Integer)::DetectionStatus = state.individuals[person_id].detected
isdetected(state::SimState, person_id::Integer)::Bool = (Detected == detected(state, person_id))

subjecthealth(state::SimState, event::Event)::HealthState = health(state, subject(event))
subjectfreedom(state::SimState, event::Event)::FreedomState = freedom(state, subject(event))

sourcehealth(state::SimState, event::Event)::HealthState = health(state, source(event))
sourcefreedom(state::SimState, event::Event)::FreedomState = freedom(state, source(event))

strainof(state::SimState, person_id::Integer) = strainof(state.forest, person_id)

#forwardinfections(state::SimState, person_id::Integer) = inclusive(state.infections, searchequalrange(state.infections, person_id)...) |> values
forwardinfections(state::SimState, person_id::Integer) = forwardinfections(state.forest, person_id)
backwardinfection(state::SimState, person_id::Integer) = backwardinfection(state.forest, person_id)


function sethealth!(state::SimState, person_id::Integer, new_health::HealthState)
  orig = state.individuals[person_id]
  @assert orig.health <= new_health
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

function quarantine_advance!(state::SimState, person_id::Integer, adv_val::Integer)
  orig = state.individuals[person_id]
  if adv_val >= 0
    state.individuals[person_id] = @set orig.quarantine_level = orig.quarantine_level+adv_val
  else
    state.individuals[person_id] = @set orig.quarantine_level = orig.quarantine_level-(-adv_val)
  end
  nothing
end

function quarantine_cancel!(state::SimState, person_id::Integer)
  orig = state.individuals[person_id]
  state.individuals[person_id] = @set orig.quarantine_level = 0
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

  event = Event(Val(OutsideTransmissionEvent), 0.0, subject(event))
  push!(state.queue, event)

end

function outsidefeed!(state::SimState, num_initial_infections::Integer, strain::StrainKind=ChineseStrain, infection_time::Real=0.0)
  N = length(state.individuals)
  individuals = 1:N

  for _ in 1:num_initial_infections
    person_id = sample(state.rng, individuals)

    event = Event(Val(OutsideInfectionEvent), infection_time, person_id, strain)
    push!(state.queue, event)
  end

end


saveparams(dict, state::SimState, prefix::AbstractString="") = saveparams(dict, state.forest, prefix)