using Random

struct AliasSampler{IndexType <: Integer, ProbabilityParamType <: Real}
    # Some alias_indices will be left uninitialized after the
    # constructor finishes. This is not a bug!
    alias_indices::Vector{IndexType}
    nonalias_probs::Vector{ProbabilityParamType}
end

function AliasSampler(IndexType::Type,
                weights::AbstractVector{ProbabilityParamType}
            ) where{ProbabilityParamType <: Real}
    @assert IndexType <: Integer

    n = length(weights)
    s = sum(weights)
    scf = n/s
    sc_probabilities = (weights*scf)::Vector{ProbabilityParamType}

    smalls = Vector{IndexType}(undef, n)
    smalls_idx = 0
    larges = Vector{IndexType}(undef, n)
    larges_idx = 0

    for idx in 1:n
        if sc_probabilities[idx] >= 1.0
            larges_idx += 1
            larges[larges_idx] = idx
        else
            smalls_idx += 1
            smalls[smalls_idx] = idx
        end
    end

    aliases = Vector{IndexType}(undef, n)

    while larges_idx > 0 && smalls_idx > 0
        sm_p_idx = smalls[smalls_idx]
        lg_p_idx = larges[larges_idx]
        aliases[sm_p_idx] = lg_p_idx

        # This is slightly better numerically than the more obvious: sc_probabilities[lg_p_idx] += sc_probabilities[sm_p_idx] - 1.0
        sc_probabilities[lg_p_idx] -= 1.0
        sc_probabilities[lg_p_idx] += sc_probabilities[sm_p_idx]

        if sc_probabilities[lg_p_idx] < 1.0
            smalls[smalls_idx] = lg_p_idx
            larges_idx -= 1
        else
            smalls_idx -= 1
        end
    end

    while larges_idx > 0
        sc_probabilities[larges[larges_idx]] = 2.0
        larges_idx -= 1
    end

    while smalls_idx > 0
        sc_probabilities[smalls[smalls_idx]] = 2.0
        smalls_idx -= 1
    end

    return AliasSampler(aliases, sc_probabilities)
end

function asample(alias_sampler::AliasSampler{IndexType, ProbabilityParamType}, rng=Random.GLOBAL_RNG)::IndexType where{IndexType <: Integer, ProbabilityParamType <: Real}

    idx = rand(rng, 1:length(alias_sampler.alias_indices))
    if alias_sampler.nonalias_probs[idx] >= 1.0
        return idx
    end
    if rand(rng) < alias_sampler.nonalias_probs[idx]
        return idx
    else
        return alias_sampler.alias_indices[idx]
    end
end

# generalized setup:

function setup_alias_sampler!(
  weights::AbstractVector{T} where T <: Real,
  acceptances::AbstractVector{T} where T <: Real,
  aliases::AbstractVector{Idx},
  smalls::AbstractVector{Idx},
  larges::AbstractVector{Idx},
  wsum::Real = sum(weights)
  ) where Idx <: Integer

  n = length(weights)
  length(acceptances) == length(aliases) == length(smalls) == length(larges) == n || throw(DimensionMismatch("Inconsistent array lengths."))

  scf = n/wsum
  for i = 1:n
    @inbounds acceptances[i] = weights[i] * scf
  end

  smalls_idx = 0
  larges_idx = 0

for idx in 1:n
    if acceptances[idx] >= 1.0
      larges_idx += 1
      larges[larges_idx] = idx
    else
      smalls_idx += 1
      smalls[smalls_idx] = idx
    end
  end

  while larges_idx > 0 && smalls_idx > 0
    sm_p_idx = smalls[smalls_idx]
    lg_p_idx = larges[larges_idx]
    aliases[sm_p_idx] = lg_p_idx

    # This is slightly better numerically than the more obvious: sc_probabilities[lg_p_idx] += sc_probabilities[sm_p_idx] - 1.0
    acceptances[lg_p_idx] -= 1.0
    acceptances[lg_p_idx] += acceptances[sm_p_idx]

    if acceptances[lg_p_idx] < 1.0
      smalls[smalls_idx] = lg_p_idx
      larges_idx -= 1
    else
      smalls_idx -= 1
    end
  end

  while larges_idx > 0
    acceptances[larges[larges_idx]] = 2.0
    larges_idx -= 1
  end

  while smalls_idx > 0
    acceptances[smalls[smalls_idx]] = 2.0
    smalls_idx -= 1
  end
  nothing
end
