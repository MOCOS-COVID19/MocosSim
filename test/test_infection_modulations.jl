#import MocosSim: StrainKind, ContactKind, ConstantKernelContact, NullStrain,
#                 TransmissionEvent, Event, HouseholdContact
using MocosSim: tanh_modulation

@testset "InfectionModulations" begin
  @testset "ReferenceFunction" begin
    @test tanh_modulation(-Inf, 5, 1, -5, 7) == -5
    @test tanh_modulation(Inf, 5, 1, -5, 7) == 7
    @test tanh_modulation(5, 5, 1, -5, 7) == (-5+7)/2

    @test tanh_modulation(-Inf, 5, 1, 3, 3) == 3
    @test tanh_modulation(0, 5, 1, 3, 3) == 3
    @test tanh_modulation(+Inf, 5, 1, 3, 3) == 3

    @testset "Randomized testing" begin
      rng = Random.MersenneTwister(13)
      for i in 1:100
        loc = rand(rng) * 200 - 100
        scale = rand(rng) * 10
        lhs = rand(rng) * 20 - 10
        rhs = rand(rng) * 20 - 10
        @test tanh_modulation(-Inf, loc, scale, lhs, rhs) ≈ lhs
        @test tanh_modulation(+Inf, loc, scale, lhs, rhs) ≈ rhs
        @test tanh_modulation(loc, loc, scale, lhs, rhs) ≈ (lhs + rhs)/2
      end
    end
  end

  struct TanhMockState <: MocosSim.AbstractSimState
    rng::MersenneTwister
  end

  MocosSim.time(::TanhMockState) = 10
  MocosSim.numdetected(::TanhMockState) = 0
  MocosSim.numdead(::TanhMockState) = 0

  struct TanhMockParams <: MocosSim.AbstractSimParams end

@testset "TanhModulation" begin
  modulation = MocosSim.TanhModulation(
      weight_detected=0,
      weight_deaths=0,
      weight_days=1,
      loc=10.0,
      scale=1.0,
      limit_value=0.0)

  mock_state = TanhMockState(MersenneTwister(13))
  strain = MocosSim.ChineseStrain

  hits = 0
  num_samples = 100
  for i in 1:num_samples
    household_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.HouseholdContact, strain)
    constant_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.ConstantKernelContact, strain)

    @test modulation(mock_state, TanhMockParams(), household_infection) == true

    hits += modulation(mock_state, TanhMockParams(), constant_infection)
  end

  @test hits == num_samples/2

end

@testset "IncreasingTanhModulation" begin
  modulation = MocosSim.IncreasingTanhModulation(
      weight_detected=0,
      weight_deaths=0,
      weight_days=1,
      loc=10.0,
      scale=1.0,
      initial_value=0.0)

  mock_state = TanhMockState(MersenneTwister(13))
  strain = MocosSim.ChineseStrain

  hits = 0
  num_samples = 100
  for i in 1:num_samples
    household_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.HouseholdContact, strain)
    constant_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.ConstantKernelContact, strain)

    @test modulation(mock_state, TanhMockParams(), household_infection) == true

    hits += modulation(mock_state, TanhMockParams(), constant_infection)
  end

  @test hits == num_samples/2

end

  #state = MocosSim.SimState(3)
  #params::DummyParams = DummyParams()
  #strain::StrainKind = NullStrain
  #contact_kind::ContactKind = ConstantKernelContact
  #event = Event(Val(TransmissionEvent), 0.0, 0, 0, contact_kind, strain)
  #@testset "TanhModulation" begin
  #  weight_detected::Float64 = 0
  #  weight_deaths::Float64 = 0
  #  weight_days::Float64 = 1
  #  loc::Float64 = 10
  #  scale::Float64 = 1
  #  limit_value::Float64 = 0.0
  #  rng = MersenneTwister(13)
#
#    modulation = MocosSim.TanhModulation(weight_detected, weight_deaths, weight_days, loc, scale,
#                                         limit_value)
#    @testset "Limit value 0.0 results in rejecting all events" begin
#      state.time = 100
#      @test modulation(state, params, event) == false
#    end
#    @testset "Household event is not modified" begin
#      household_event = Event(Val(TransmissionEvent), 0.0, 0, 0, HouseholdContact, strain)
#      @test modulation(state, params, household_event) == true
#    end
#  end
end
