using MocosSim
using CodecZlib
using FileIO
using JLD2
using DataFrames
using FixedPointNumbers
using Random
using Profile


individuals_df = load( (@__DIR__) * "/wroclaw_v4.jld2", "individuals_df")
rng = Random.MersenneTwister(0);

@info "creating state"
@time state = MocosSim.SimState(nrow(individuals_df))

@info "loading data"
@time params = MocosSim.load_params(
  population=individuals_df,
  mild_detection_prob=0.4,
  backward_tracing_prob=0.2,
  forward_tracing_prob=0.2,
  constant_kernel_param=0.7425,
  household_kernel_param=0.07,
  british_strain_multiplier=1.7,

  age_coupling_thresholds=[0, 18, 40],
  age_coupling_weights=[1.1 1.2 1.3;
                        2.1 2.2 2.3;
                        3.1 3.2 3.3],
  age_coupling_use_genders=false,
  age_coupling_param=1.0,

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
  testing_time=3.0,
  hospitalization_time_ratio = 0.5,
  hospitalization_multiplier=1.5,
  death_multiplier = 0.5
);

mutable struct Callback
  max_time::Float64
  num_events::Int
end

function (cb::Callback)(event::MocosSim.Event, ::MocosSim.SimState, ::MocosSim.SimParams)::Bool
  cb.num_events += 1
  MocosSim.time(event) < cb.max_time
end

@info "initializing state"
@time MocosSim.reset!(state, MersenneTwister(0))
@time MocosSim.InstantOutsideCases(;num_infections=100)(state, params)

immunity_params = MocosSim.make_immunity_table(state, 0.5) 
MocosSim.immunize!(state,immunity_params; enqueue = true)

@info "warm-up"
cb = Callback(50, 0)
@time MocosSim.simulate!(state, params, cb)
@info "events executed = $(cb.num_events), simtime=$(state.time)"

@info "measuring allocations"
cb = Callback(100, 0)
Profile.clear_malloc_data()
@time MocosSim.simulate!(state, params, cb)
@info "events executed = $(cb.num_events)"
