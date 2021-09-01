import MocosSim: StrainKind, ContactKind, ConstantKernelContact, NullStrain,
                 TransmissionEvent, Event, HouseholdContact

mutable struct DummyParams <: MocosSim.AbstractSimParams
  DummyParams() = new()
end

@testset "InfectionModulations" begin
  state = MocosSim.SimState(3)
  params::DummyParams = DummyParams()
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
    @testset "Limit value 0.0 results in rejecting all events" begin
      state.time = 100
      @test modulation(state, params, event) == false
    end
    @testset "Household event is not modified" begin
      household_event = Event(Val(TransmissionEvent), 0.0, 0, 0, HouseholdContact, strain)
      @test modulation(state, params, household_event) == true
    end
  end
end
