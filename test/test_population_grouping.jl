@testset "PopulationGrouping" begin
  @testset "ExampleGrouping" begin
    num_groups = 100
    num_individuals = 100
    rng = MersenneTwister(13)
    group_ids = rand(rng, 1:num_groups, num_individuals)
    @test minimum(group_ids) > 0
    @test maximum(group_ids) <= num_groups

    grouping = MocosSim.PopulationGrouping(group_ids, num_groups)

    @test MocosSim.numgroups(grouping) == num_groups

    @testset "Every group has the right people" begin
      for group_id in 1:num_groups
        group = MocosSim.getgroup(grouping, group_id)

        for person_id in group
          @test group_ids[person_id] == group_id
        end
      end
    end

    @testset "Every person is in the right group" begin
      for person_id in 1:num_individuals
        group_id = group_ids[person_id]
        group = MocosSim.getgroup(grouping, group_id)
        @test person_id in group
      end
    end

    @testset "Group sizes" begin
      sizes = MocosSim.groupsizes(grouping)
      @test length(sizes) == num_groups
      for i in 1:MocosSim.numgroups(grouping)
        @test length(MocosSim.getgroup(grouping, i)) == sizes[i]
        @test MocosSim.groupsize(grouping, i) == sizes[i]
      end
    end
  end
end