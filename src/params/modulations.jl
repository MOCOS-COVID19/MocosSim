struct TanhModulation <: InfectionModulation
  weight_detected::Float64
  weight_deaths::Float64
  weight_days::Float64
  loc::Float64
  scale::Float64
  limit_value::Float64

  TanhModulation(;weight_detected::Real=0, weight_deaths::Real=0, weight_days::Real=0, loc::Real=0, scale::Real=0, limit_value::Real=0) =
    0 <= limit_value <= 1 ? new(weight_detected, weight_deaths, weight_days, loc, scale, limit_value) : error("limit_value must be from 0 to 1, got $limit_value")
end

function tanh_modulation(x::Real, loc::Real, scale::Real, lhs_limit::Real, rhs_limit::Real)
  t_11 = tanh((x - loc)/scale)
  t_01 = (t_11+1)/2
  t_01 * (rhs_limit-lhs_limit) + lhs_limit
end

function evalmodulation(f::TanhModulation, state::AbstractSimState, ::AbstractSimParams)::Float64
  num_days = time(state)
  num_detected = numdetected(state)
  num_deaths = numdead(state)
  fear = num_detected * f.weight_detected + num_deaths * f.weight_deaths + num_days * f.weight_days
  tanh_modulation(fear, f.loc, f.scale, 1.0, f.limit_value)
end

struct IncreasingTanhModulation <: InfectionModulation
  weight_detected::Float64
  weight_deaths::Float64
  weight_days::Float64
  loc::Float64
  scale::Float64
  initial_value::Float64

  IncreasingTanhModulation(;weight_detected::Real=0, weight_deaths::Real=0, weight_days::Real=0, loc::Real=0, scale::Real=0, initial_value::Real=0) =
    0 <= initial_value <= 1 ? new(weight_detected, weight_deaths, weight_days, loc, scale, initial_value) : error("initial_value must be from 0 to 1, got $initial_value")
end

function evalmodulation(f::IncreasingTanhModulation, state::AbstractSimState, ::AbstractSimParams)::Float64
  num_days = time(state)
  num_detected = numdetected(state)
  num_deaths = numdead(state)

  fear = num_detected * f.weight_detected + num_deaths * f.weight_deaths + num_days * f.weight_days
  tanh_modulation(fear, f.loc, f.scale, f.initial_value, 1.0)
end

evalmodulation(modulation::Nothing, state::AbstractSimState, params::AbstractSimParams)::Float64 = 1.0

function infectionsuccess(modulation::InfectionModulation, state::AbstractSimState, params::AbstractSimParams, event::Event)::Bool
  @assert kind(event) == TransmissionEvent

  ck = contactkind(event)
  if ConstantKernelContact !== ck && AgeCouplingContact !== ck
    return true # do not affect other types of contact than "outer" ones
  end

  rand(state.rng) < evalmodulation(modulation, state, params)
end

function infectionsuccess(state::AbstractSimState, params::AbstractSimParams, event::Event)::Bool
  if isnothing(params.infection_modulation)
    return true
  end

  @assert kind(event) == TransmissionEvent

  ck = contactkind(event)
  if ConstantKernelContact !== ck && AgeCouplingContact !== ck
    return true # do not affect other types of contact than "outer" ones
  end

  rand(state.rng) < evalmodulation(params.infection_modulation, state, params)
#  infectionsuccess(params.infection_modulation, state, params, event)
end

# This all to avoid using @eval and others
const modulations = Dict{String, Type{T} where T}(
    "TanhModulation" => TanhModulation,
    "IncreasingTanhModulation" => IncreasingTanhModulation
)

make_infection_modulation(::Nothing) = nothing
make_infection_modulation(name::AbstractString; args...) = modulations[name](;args...)
