@testset "HouseholdGrouping" begin
  @testset "groupptrs" begin
    household_flags = [1, 1, 2, 3, 5, 5, 5, 6, 8, 9, 10]
    begs, ends = MocosSim.groupptrs(household_flags)

    num_individuals = length(household_flags)

    for i in 1:num_individuals
      household = begs[i]:ends[i]
      for j in 1:num_individuals
        if household_flags[i] == household_flags[j]
          @test j ∈ household
        else
          @test j ∉ household
        end
      end
    end
  end
end