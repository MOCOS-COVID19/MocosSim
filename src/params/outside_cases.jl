Base.@kwdef struct ParabolicOutsideCases <: AbstractOutsideCases
  days::Float64 = 30
  peak::Float64 = 20
  height::Float64 = 0.8
  minimum::Float64 = 0.02
  frequency::TimePoint = 0.05 |> MocosSim.TimePoint
  strain::StrainKind = ChineseStrain
  time_limit::TimePoint = typemax(TimePoint)
end
  
function (f::ParabolicOutsideCases)(state::AbstractSimState, ::AbstractSimParams)
  days = (f.days/2)^2
  N = length(state.individuals)
  individuals = 1:N

  for infection_time in 0.0:f.frequency:time_limit
    if rand(state.rng) < max((Float64(-abs(infection_time-f.peak)^2)*f.height)/days + f.height, f.minimum)
      person_id = sample(state.rng, individuals)
      event = Event(Val(OutsideInfectionEvent), infection_time, person_id, strain)
      push!(state.queue, event)
    end
  end

end

Base.@kwdef struct InstantOusideCases <: AbstractOutsideCases
  num_infections::Float64 = 50
  import_time::TimePoint = 0.0 |> MocosSim.TimePoint
  strain::StrainKind = ChineseStrain
end
  
function (f::InstantOusideCases)(state::AbstractSimState, ::AbstractSimParams)
  N = length(state.individuals)
  individuals = 1:N

  for _ in 1:f.num_infections
    person_id = sample(state.rng, individuals)
    event = Event(Val(OutsideInfectionEvent), f.import_time, person_id, strain)
    push!(state.queue, event)
  end
end

# This all to avoid using @eval and others
const imported_cases = Dict{String, Type{T} where T}(
    "ParabolicOutsideCases" => ParabolicOutsideCases,
    "InstantOusideCases"   => InstantOusideCases
)

make_imported_cases(name::AbstractString; args...) = imported_cases[name](;args...)
