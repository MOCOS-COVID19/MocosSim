module MocosSim

using Distributions: minimum
using DataStructures: minimum
using LinearAlgebra: AbstractMatrix
__precompile__(true)

using DataFrames
using DataStructures
using Distributions
using FixedPointNumbers
using LinearAlgebra
using Random
using SaferIntegers
using Setfield
using StaticArrays

import Base.show
import Distributions.sample

const TimePoint = Fixed{Int32, 16}
const TimeDiff = Fixed{Int32, 16}
const TimeDay = Int16
const PersonIdx=UInt32

include("utils.jl")

include("enums.jl")
include("event.jl")

include("eventqueue.jl")
include("robin_forest.jl")
const InfectionForest = RobinForest
include("simstate.jl")

include("simparams.jl")

include("event_execution.jl")
include("infection_kernels.jl")

export simulate!
export Event

function simulate!(state::SimState, params::SimParams)
  while true
    if isempty(state.queue)
      break
    end
    event = pop!(state.queue)
    execute!(state, params, event)
  end
  nothing
end

function simulate!(state::SimState, params::SimParams, callback)
  while true
    if isempty(state.queue)
      break
    end

    event = pop!(state.queue)

    executed::Bool = execute!(state, params, event)

    if executed
       should_continue = callback(event, state, params)
       if !should_continue
         break
       end
    end
  end
  nothing
end

precompile(simulate!, (SimState, SimParams))

end #module
