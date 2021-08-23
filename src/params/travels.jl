Base.@kwdef struct Holiday <: InfectionTravels
    days::Float64 = 0
    peak::Float64 = 0
    height::Float64 = 0
    minimum::Float64 = 0
  end
  
  function (f::Holiday)(state::SimState, params::SimParams, event::Event, infection_time::Real=0.0)
    @assert kind(event) == OutsideTransmissionEvent
    days = (f.days/2)^2
    rand(state.rng) < max((Float64(-abs(event.time-f.peak)^2)*f.height)/days + f.height, f.minimum)
  end
  
  # This all to avoid using @eval and others
  const travels = Dict{String, Type{T} where T}(
      "Holiday" => Holiday
  )
  
  make_infection_travels(name::AbstractString; args...) = travels[name](;args...)