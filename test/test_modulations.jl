include("../src/enums.jl")
include("../src/event.jl")

@testset "InfectionModulations" begin
  state = MocosSim.SimState(3)
  params = missing
  strain::StrainKind = NullStrain
  contact_kind::ContactKind = ConstantKernelContact
  event = Event(Val(TransmissionEvent), 0.0, 0, 0, contact_kind, strain)
  @testset "TanhModulation" begin
    weight_detected::Float64 = 0
    weight_deaths::Float64 = 0
    weight_days::Float64 = 1
    loc::Float64 = 10
    scale::Float64 = 1
    limit_value::Float64 = 0.0
    rng = MersenneTwister(13)

    modulation = MocosSim.TanhModulation(weight_detected, weight_deaths, weight_days, loc, scale,
                                         limit_value)
    state.time = 100
    @test modulation(state, params, event) == false

    household_event = Event(Val(TransmissionEvent), 0.0, 0, 0, HouseholdContact, strain)
    @test modulation(state, params, household_event) == true
  end
end
