using MocosSim
using FileIO
using JLD2
using DataFrames
using FixedPointNumbers
using Random
using ProgressMeter
using PyPlot
using StatsBase

using Profile


individuals_df = load("../dane/poland_v3.jld2", "individuals_df");
rng = Random.MersenneTwister(0);

state = MocosSim.SimState(nrow(individuals_df))

params = MocosSim.load_params(
  population=individuals_df,
  mild_detection_prob=0.4,
  backward_tracing_prob=0.2,
  forward_tracing_prob=0.2,
  constant_kernel_param=0.7425,
  household_kernel_param=0.07,
  british_strain_multiplier=1.7,

  infection_modulation_name="TanhModulation",
  infection_modulation_params=(
    scale=12500,
    loc=300000,
    weight_detected=1,
    weight_deaths=0,
    limit_value=0.545455
  ),

  forward_detection_delay=1.75,
  backward_detection_delay=1.75,
  testing_time=3.0
);

struct Callback
  max_time::Float64
end

function (cb::Callback)(event::MocosSim.Event, state::MocosSim.SimState, ::MocosSim.SimParams)
  MocosSim.time(event) < cb.max_time
end

MocosSim.reset!(state, MersenneTwister(0))
MocosSim.initialfeed!(state, 100)
MocosSim.simulate!(state, params, Callback(100))

Profile.clear_malloc_data()
MocosSim.simulate!(state, params, Callback(100))
