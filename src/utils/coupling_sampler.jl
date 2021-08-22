struct CouplingSampler{Tprob<:Real, PersonIdx<:Integer, GroupIdx<:Integer}
  group_ids::Vector{GroupIdx} # TODO: needed only for reverse grouping computation, forward can be computed on-the-fly
  grouping::PopulationGrouping{PersonIdx}
  matrix_sampler::MatrixAliasSampler{GroupIdx, Float64}
end

function CouplingSampler(group_ids::AbstractVector{GroupIdx}, coupling_weights::Matrix{T} where T<:Real, PersonIdx::DataType=UInt32) where GroupIdx<:Integer
  num_groups, M = size(coupling_weights)
  @assert num_groups == M

  @assert maximum(group_ids) <= num_groups <= typemax(GroupIdx)
  @assert minimum(group_ids) > 0

  @assert maximum(coupling_weights) < Inf
  @assert minimum(coupling_weights) >= 0

  grouping = PopulationGrouping(group_ids, num_groups, PersonIdx)
  matrix_sampler = MatrixAliasSampler(coupling_weights, GroupIdx)
  CouplingSampler{eltype(coupling_weights), PersonIdx, GroupIdx}(group_ids, grouping, matrix_sampler)
end

groupof(s::CouplingSampler, person_idx::Integer) = s.group_ids[person_idx]
getgroup(coupling::CouplingSampler, group_idx::Integer) = getgroup(coupling.grouping, group_idx)

function sample(rng::AbstractRNG, coupling::CouplingSampler, source_id::Integer)
  # retrieve the source group
  source_group_id = groupof(coupling, source_id)

  # alias sample the group id
  target_group_id = sample(rng, coupling.matrix_sampler, source_group_id)

  # retrieve the group
  target_group = getgroup(coupling.grouping, target_group_id)
  @assert length(target_group) > 0 "empty population group $target_group_id, grouping = $(coupling.grouping)"

  # sample from the group
  rand(rng, target_group)
end