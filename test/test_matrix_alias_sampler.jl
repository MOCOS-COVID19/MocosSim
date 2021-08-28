@testset "MatrixAliasSampler" begin
  @testset "Rates similar to empirical results" begin
    rng = MersenneTwister(13)
    N = 5
    example_rates = example_probs = rand(rng, N, N) * 10
    matrix_sampler = MocosSim.MatrixAliasSampler(example_probs)

    num_samples = 10^7
    counts = zeros(Int, N, N)

    for src in 1:N
      for _ in 1:num_samples
        tgt = MocosSim.sample(rng, matrix_sampler, src)
        counts[tgt, src] += 1
      end
    end

    empirial_probs = counts ./ num_samples
    theoretical_probs = example_rates ./ sum(example_rates, dims=1)

    @test maximum(abs.(empirial_probs .- theoretical_probs)) < 3/sqrt(num_samples)


  end
end