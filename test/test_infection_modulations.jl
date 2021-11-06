using MocosSim: tanh_modulation, infectionsuccess

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

  mutable struct TanhMockState <: MocosSim.AbstractSimState
    rng::MersenneTwister
    time::Float64
  end

  MocosSim.time(s::TanhMockState) = s.time
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

  for t in -3:0.1:3
    mock_state = TanhMockState(MersenneTwister(13), t)
    @test MocosSim.evalmodulation(modulation, mock_state, TanhMockParams()) == tanh_modulation(t, 10, 1, 1.0, 0.0)
  end

  strain = MocosSim.ChineseStrain

  hits = 0
  num_samples = 100

  mock_state = TanhMockState(MersenneTwister(13), 10.0)
  for i in 1:num_samples
    household_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.HouseholdContact, strain)
    constant_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.ConstantKernelContact, strain)

    @test infectionsuccess(modulation, mock_state, TanhMockParams(), household_infection) == true

    hits += infectionsuccess(modulation, mock_state, TanhMockParams(), constant_infection)
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

  for t in -3:0.1:3
    mock_state = TanhMockState(MersenneTwister(13), t)
    @test MocosSim.evalmodulation(modulation, mock_state, TanhMockParams()) == tanh_modulation(t, 10, 1, 0.0, 1.0)
  end

  mock_state = TanhMockState(MersenneTwister(13), 10.0)
  strain = MocosSim.ChineseStrain

  hits = 0
  num_samples = 100
  for i in 1:num_samples
    household_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.HouseholdContact, strain)
    constant_infection = Event(Val(MocosSim.TransmissionEvent), 0.0, 0, 0, MocosSim.ConstantKernelContact, strain)

    @test infectionsuccess(modulation, mock_state, TanhMockParams(), household_infection) == true

    hits += infectionsuccess( modulation, mock_state, TanhMockParams(), constant_infection)
  end

  @test hits == num_samples/2

end

end
