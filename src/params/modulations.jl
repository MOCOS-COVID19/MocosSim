Base.@kwdef struct TanhModulation <: InfectionModulation
  weight_detected::Float64 = 0
  weight_deaths::Float64 = 0
  weight_days::Float64 = 0
  loc::Float64
  scale::Float64
  limit_value::Float64
end

function (f::TanhModulation)(state::SimState, ::SimParams, event::Event)
  @assert kind(event) == TransmissionEvent

  ck = contactkind(event)
  if ConstantKernelContact !== ck && AgeCouplingContact !== ck
    return true # do not affect other types of contact than "outer" ones
  end

  num_detected = numdetected(state.stats)
  num_deaths = numdead(state.stats)
  num_days = time(state)

  fear = num_detected * f.weight_detected + num_deaths * f.weight_deaths + num_days * f.weight_days
  x = (fear - f.loc) / f.scale
  scaling = ((1 - f.limit_value) / 2)
  base = (1 - (1 - f.limit_value) / 2)
  rand(state.rng) < (-tanh(x) * scaling + base) / max(1, f.limit_value)
end

# This all to avoid using @eval and others
const modulations = Dict{String, Type{T} where T}(
    "TanhModulation" => TanhModulation
)

make_infection_modulation(name::AbstractString; args...) = modulations[name](;args...)
