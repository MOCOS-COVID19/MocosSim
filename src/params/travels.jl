Base.@kwdef struct Holiday <: InfectionTravels
    days::Float64 = 0
    peak::Float64 = 0
    height::Float64 = 0
    minimum::Float64 = 0
  end
  
  function (f::Holiday)(state::SimState, params::SimParams, event::Event, infection_time::Real=0.0)
    @assert kind(event) == OutsideTransmissionEvent
    ck = contactkind(event)
    #if ConstantKernelContact !== ck && SporadicContact !== ck && FriendshipContact !== ck
    #  return true # do not affect other types of contact than "outer" ones
    #end
    #rand(state.rng) < max((Float64(-abs(infection_time-f.peak)*abs(infection_time-f.peak)))/1500.0 + f.height,f.minimum)
    vaule = rand(state.rng)
    #probability =  max((Float64(-abs(infection_time-f.peak)*abs(infection_time-f.peak)))/1500.0 + f.height,f.minimum)
    probability =  max((Float64(-abs(event.time-f.peak)*abs(event.time-f.peak)))/1500.0 + f.height,f.minimum)
    #@info "value" vaule
    #@info "probability" probability
    vaule < probability
    
  end
  
  # This all to avoid using @eval and others
  const travels = Dict{String, Type{T} where T}(
      "Holiday" => Holiday
  )
  
  make_infection_travels(name::AbstractString; args...) = travels[name](;args...)