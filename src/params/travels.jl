Base.@kwdef struct ParabolicImports <: ImportedCases
    days::Float64 = 30
    peak::Float64 = 20
    height::Float64 = 0.8
    minimum::Float64 = 0.02
  end
  
  function (f::ParabolicImports)(state::SimState, params::SimParams, infection_time::Real=0.0)
    #@assert kind(event) == OutsideTransmissionEvent
    days = (f.days/2)^2
    rand(state.rng) < max((Float64(-abs(infection_time-f.peak)^2)*f.height)/days + f.height, f.minimum)
  end
  
  # This all to avoid using @eval and others
  const imported_cases = Dict{String, Type{T} where T}(
      "ParabolicImports" => ParabolicImports
  )
  
  make_imported_cases(name::AbstractString; args...) = imported_cases[name](;args...)
