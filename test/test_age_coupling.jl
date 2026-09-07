@testset "AgeCouplingParams" begin
  @testset "agegroup() assings to the right groups" begin

    age_thresholds = Int[0, 5, 12, 18, 30, 40, 50, 60, 70]
    reference_implementation(age) =
      age < 5 ? 1 :
      age < 12 ? 2 :
      age < 18 ? 3 :
      age < 30 ? 4 :
      age < 40 ? 5 :
      age < 50 ? 6 :
      age < 60 ? 7 :
      age < 70 ? 8 :
      9

      for age in 0:100
        @test MocosSim.agegroup(age_thresholds, age) == reference_implementation(age)
      end

  end
end


@testset "gender-aware AgeCouplingParams" begin
  ages = [4, 8, 4, 8]
  genders = Bool[false, false, true, true]
  thresholds = [0, 5]
  weights = [1.0 0 0 0; 0 1.0 0 0; 0 0 1.0 0; 0 0 0 1.0]
  params = MocosSim.AgeCouplingParams(ages, genders, thresholds, weights, nothing)
  @test length(params.source_weighting) == 4
  @test params.coupling.group_ids == MocosSim.GroupIdx[1, 3, 2, 4]
  @test_throws DimensionMismatch MocosSim.AgeCouplingParams(
    ages, genders, thresholds, Matrix{Float64}(undef, 2, 2), nothing)

  params_without_genders = MocosSim.AgeCouplingParams(
    ages, nothing, thresholds, [1.0 0; 0 1.0], nothing)
  @test params_without_genders.coupling.group_ids == MocosSim.GroupIdx[1, 2, 1, 2]
end
