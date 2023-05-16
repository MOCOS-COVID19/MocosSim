const GroupIdx = UInt8

struct AgeCouplingParams
  coupling::CouplingSampler{Float64, PersonIdx, GroupIdx}
  source_weighting::Vector{Float64}
end

agegroup(age_thresholds::AbstractArray{T} where T<:Real, age::Real) = GroupIdx(searchsortedlast(age_thresholds, age))

function AgeCouplingParams(
  ages::AbstractVector{T} where T<:Real,
  genders::Union{Nothing, AbstractVector{Bool}},
  age_thresholds::AbstractVector{T} where T<:Real,
  coupling_weights::Matrix{Float64},
  contact_normalization::Union{Nothing, Real})

  @assert age_thresholds[1] == 0
  num_groups = length(age_thresholds)
  @assert (num_groups,num_groups) == size(coupling_weights)

  group_ids = agegroup.((age_thresholds,), ages)
  @assert minimum(group_ids) > 0
  if genders !== nothing
    @assert length(ages) == length(genders)
    group_ids .= group_ids .*2 .+ genders .- 1
  end
  @assert maximum(group_ids) <= num_groups <= typemax(GroupIdx) "minmax: $(extrema(group_ids)), num=$num_groups, typemax=$(typemax(GroupIdx))"

  source_weighting = Vector{Float64}(undef, size(coupling_weights, 2))
  sum!(source_weighting', coupling_weights)

  coupling_sampler = CouplingSampler(group_ids, coupling_weights, PersonIdx)

  if nothing !== contact_normalization # normalize such that expected number of contacts is equal to the given value
    group_sizes = groupsizes(coupling_sampler.grouping)
    mean_contacts = source_weighting' * group_sizes / sum(group_sizes)
    source_weighting .*= contact_normalization / mean_contacts
  end
  AgeCouplingParams(
    coupling_sampler,
    source_weighting
  )
end

sourceweightof(p::AgeCouplingParams, source_id::Integer) = p.source_weighting[groupof(p.coupling, source_id)]