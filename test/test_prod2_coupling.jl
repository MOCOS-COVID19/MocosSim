@testset "Prod2Coupling" begin
  @testset "jointgrouping is indexing in the right order" begin
    num_groups1 = 7
    num_groups2 = 13

    joint_ids = Int[]
    for group1_id in 1:num_groups1
      for group2_id in 1:num_groups2
        joint_group_id = MocosSim.jointgrouping(group2_id, num_groups2, group1_id)
        push!(joint_ids, joint_group_id)
      end
    end
    @test issorted(joint_ids)
    @test length(unique(joint_ids)) == num_groups1 * num_groups2
  end

  @testset "sampling with right probability" begin
    function empiricalprobs(rng::AbstractRNG, coupling::MocosSim.Prod2CouplingSampler,
        source_g1::Integer, num_g1::Integer, source_g2::Integer, num_g2::Integer,
        group1_ids::AbstractVector{T} where T<:Integer, group2_ids::AbstractVector{T} where T<:Integer;
        N::Int=10^6)
      num_individuals = length(group1_ids)
      sampled_people = UInt32[]
      for _ in 1:N
          person_id = MocosSim.sample(rng, coupling, source_g1, source_g2)
          push!(sampled_people, person_id)
      end
      hist = fit(Histogram, sampled_people, 0:num_individuals, closed=:right)

      df = DataFrame(
        g1 = group1_ids,
        g2 = group2_ids,
        g = MocosSim.jointgrouping.(group2_ids, num_g2, group1_ids),
          results = hist.weights / N
      )
     df = combine(groupby(df, :g), :results => sum => :prob)
     df = leftjoin(DataFrame(g=1:(num_g1*num_g2)), df, on=:g)
     df = sort(df, :g)
     replace!(df.prob, missing => 0.0) .|> Float64
    end

    @testset "all coupling equal" begin
      rng = MersenneTwister(13)

      num_groups1 = 3
      num_groups2 = 5

      group1_ids = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3]
      group2_ids = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 1, 2, 3, 5, 4]

      weights1 = ones(num_groups1, num_groups1)
      weights2 = ones(num_groups2, num_groups2)

      coupling = MocosSim.Prod2CouplingSampler(group1_ids, weights1, group2_ids, weights2)

      prod_weights = kron(weights1, weights2);
      prod_weights .*= MocosSim.groupsizes(coupling.grouping)
      prod_weights ./= sum(prod_weights, dims=1);

      for g1 in 1:num_groups1
        for g2 in 1:num_groups2
          empirical = empiricalprobs(rng, coupling, g1, num_groups1, g2, num_groups2, group1_ids, group2_ids)
          theoretical = prod_weights[:, MocosSim.jointgrouping.(g2, num_groups2, g1)]
          @test norm(empirical .- theoretical, Inf) < 0.001

        end
      end
    end
  end
end