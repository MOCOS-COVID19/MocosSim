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