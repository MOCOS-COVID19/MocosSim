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

struct TwoTanhModulations <: InfectionModulation
  weight_detected::Float64
  weight_deaths::Float64
  weight_days::Float64
  loc1::Float64
  scale1::Float64
  limit_value1::Float64
  loc2::Float64
  scale2::Float64
  limit_value2::Float64

  TwoTanhModulations(;weight_detected::Real=0, weight_deaths::Real=0, weight_days::Real=0, loc1::Real=0, scale1::Real=0, limit_value1::Real=0, loc2::Real=0, scale2::Real=0, limit_value2::Real=0) =
    0 <= limit_value2 <= limit_value1 <= 1 && loc1 < loc2 ? new(weight_detected, weight_deaths, weight_days, loc1, scale1, limit_value1, loc2, scale2, limit_value2) : error("initial_value1 and initial_value2 must be from 0 to 1 and initial_value2 must be not larger than initial_value1, got ($initial_value1, $initial_value2), also loc2 must be after loc1, got ($loc1, $loc2)")
end

function evalmodulation(f::TwoTanhModulations, state::AbstractSimState, ::AbstractSimParams)::Float64
  num_days = time(state)
  num_detected = numdetected(state)
  num_deaths = numdead(state)
  fear = num_detected * f.weight_detected + num_deaths * f.weight_deaths + num_days * f.weight_days
  modulation1 = tanh_modulation(fear, f.loc1, f.scale1, 1.0, f.limit_value1)
  tanh_modulation(fear, f.loc2, f.scale2, modulation1, f.limit_value2)
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

struct IncreasingTwoTanhModulations <: InfectionModulation
  weight_detected::Float64
  weight_deaths::Float64
  weight_days::Float64
  loc1::Float64
  scale1::Float64
  initial_value1::Float64
  loc2::Float64
  scale2::Float64
  initial_value2::Float64

  IncreasingTwoTanhModulations(;weight_detected::Real=0, weight_deaths::Real=0, weight_days::Real=0, loc1::Real=0, scale1::Real=0, initial_value1::Real=0, loc2::Real=0, scale2::Real=0, initial_value2::Real=0) =
    0 <= initial_value1 <= initial_value2 <= 1 && loc1 < loc2 ? new(weight_detected, weight_deaths, weight_days, loc1, scale1, initial_value1, loc2, scale2, initial_value2) : error("initial_value1 and initial_value2 must be from 0 to 1 and initial_value2 must be not smaller than initial_value1, got ($initial_value1, $initial_value2), also loc2 must be after loc1, got ($loc1, $loc2)")
end

function evalmodulation(f::IncreasingTwoTanhModulations, state::AbstractSimState, ::AbstractSimParams)::Float64
  num_days = time(state)
  num_detected = numdetected(state)
  num_deaths = numdead(state)

  fear = num_detected * f.weight_detected + num_deaths * f.weight_deaths + num_days * f.weight_days
  modulation1 = tanh_modulation(fear, f.loc1, f.scale1, f.initial_value1, f.initial_value2)
  tanh_modulation(fear, f.loc2, f.scale2, modulation1, 1.0)
end

struct IntervalsModulations <: InfectionModulation
  weight_detected::Float64
  weight_deaths::Float64
  weight_days::Float64
  interval_values::Vector{Float64}
  interval_times::Vector{TimePoint}

  IntervalsModulations(;weight_detected::Real=0, weight_deaths::Real=0, weight_days::Real=0, interval_values::Vector{Real}=Float64[1.0], interval_times::Vector{Real}=TimePoint[]) =
    length(interval_values) == length(interval_times) + 1 ? new(weight_detected, weight_deaths, weight_days, interval_values, interval_times) : error("length of interval_values vector must be one longer than interval_times")
end

function intervals_modulation(x::Real, interval_values::Vector{Real}, interval_times::Vector{Real})
  idx = searchsortedlast(insert!(interval_times,1,0.0), x)
  interval_values[idx]
end

function evalmodulation(f::IntervalsModulations, state::AbstractSimState, ::AbstractSimParams)::Float64
  num_days = time(state)
  num_detected = numdetected(state)
  num_deaths = numdead(state)
  fear = num_detected * f.weight_detected + num_deaths * f.weight_deaths + num_days * f.weight_days
  intervals_modulation(fear, f.interval_values, f.interval_times)
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
    "TwoTanhModulations" => TwoTanhModulations,
    "IncreasingTanhModulation" => IncreasingTanhModulation,
    "IncreasingTwoTanhModulations" => IncreasingTwoTanhModulations,
    "IntervalsModulations" => IntervalsModulations
)

make_infection_modulation(::Nothing) = nothing
make_infection_modulation(name::AbstractString; args...) = modulations[name](;args...)
