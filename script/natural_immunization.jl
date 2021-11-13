module NaturalImmunization
  push!(LOAD_PATH, "__DIR__"*"../MocosSim")
  using MocosSim
  using Random
  using JLD2
  using FileIO
  using DataFrames
  using CodecZlib
  using FixedPointNumbers
  using ArgParse

  function run(ARGS)
    parser = ArgParseSettings()
    @add_arg_table! parser begin
      "POPULATION"
        help = "path to the popultion file"
        required = true
        arg_type = String
      "IMMUNIZATION"
        help = "where to save the immunization"
        required = true
        arg_type = String
    end
    args = parse_args(ARGS, parser)


    individuals_df = load(args["POPULATION"], "individuals_df")

    N = nrow(individuals_df)
    state = MocosSim.SimState(N)

    params = MocosSim.load_params(
      population=individuals_df,
      mild_detection_prob=0.3,
      backward_tracing_prob=0.2,
      forward_tracing_prob=0.2,
      constant_kernel_param=1.35,
      household_kernel_param=0.07
    );

    MocosSim.reset!(state, MersenneTwister(0))
    MocosSim.initialfeed!(state, 100)
    MocosSim.simulate!(state, params)

    times = MocosSim.time.(state.forest.inedges)
    times[MocosSim.kind.(state.forest.inedges) .== MocosSim.InvalidEvent] .= typemax(Fixed{Int32,16})
    ordering = times |> sortperm |> invperm

    save(args["IMMUNIZATION"], "ordering", ordering .|> UInt32, compress=true)
  end

  precompile(run, (Vector{String},) )
end

NaturalImmunization.run(ARGS)