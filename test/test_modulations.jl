
@testset "InfectionModulations" begin
  struct MockSimParams <: MocosSim.AbstractSimParams end
  struct MockSimState <: MocosSim.AbstractSimState
    rng
  end

  MocosSim.time(::MockSimState) = 100
  MocosSim.numdead(::MockSimState) = 0
  MocosSim.numdetected(::MockSimState) = 0

  state = MockSimState(MersenneTwister(13))
  params = MockSimParams()

  strain::MocosSim.StrainKind = MocosSim.NullStrain
  contact_kind::MocosSim.ContactKind = MocosSim.ConstantKernelContact
  event = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, contact_kind, strain)
  @testset "TanhModulation" begin
    modulation = MocosSim.TanhModulation(
      weight_detected=0,
      weight_deaths=0,
      weight_days=1,
      loc=10.0,
      scale=1.0,
      limit_value=0.0)

    @test modulation(state, params, event) == false

    household_event = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.HouseholdContact, strain)
    @test modulation(state, params, household_event) == true
  end
end
