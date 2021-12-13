# compact end efficient sampler for 2D distribution

struct MatrixAliasSampler{IndexType<:Integer, ProbType<:Real}
  acceptances_mat::Matrix{ProbType}
  aliases_mat::Matrix{IndexType}
  # each column of both matrices represent an alias sampler
end

function MatrixAliasSampler(probs::AbstractMatrix{T}, IdxType::Type=Int) where T<:Real
  n = size(probs, 1)
  @assert size(probs, 2) == n
  aliases_mat = zeros(IdxType, n, n)
  acceptances_mat = zeros(T, n, n)

  smalls = Vector{IdxType}(undef, n)
  larges = Vector{IdxType}(undef, n)

  for i in 1:n
    weights = view(probs, :, i)
    aliases = view(aliases_mat, :, i)
    acceptances = view(acceptances_mat, :, i)

    setup_alias_sampler!(weights, acceptances, aliases, smalls, larges, sum(weights))
  end
  MatrixAliasSampler{IdxType, T}(acceptances_mat, aliases_mat)
end

function sample(rng::AbstractRNG, m::MatrixAliasSampler, source::Integer)
  N, M = m.aliases_mat |> size

  @assert N==M
  @assert N == size(m.acceptances_mat,1) == size(m.acceptances_mat, 2)
  acceptances = @view m.acceptances_mat[:, source]
  aliases = @view m.aliases_mat[:, source]

  Idx = eltype(aliases)
  idx = rand(rng, Idx(1):Idx(N))
  if acceptances[idx] >= 1.0
    return idx
  end
  rand(rng) < acceptances[idx] ? idx : aliases[idx]
end