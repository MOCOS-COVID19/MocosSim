struct Prod2CouplingSampler{GroupIdx1<:Integer, GroupIdx2<:Integer, PersonIdx<:Integer}
  group_ids1::Vector{GroupIdx1}
  matrix_sampler1::MatrixAliasSampler{GroupIdx1, Float64}

  group_ids2::Vector{GroupIdx2}
  weights_matrix2::Matrix{Float64}

  grouping::PopulationGrouping{PersonIdx}
end

jointgrouping(g1::Integer, num_g1::Integer, g2::Integer) = 1 + (g1-1) + (g2-1) * num_g1

function Prod2CouplingSampler(
    group_ids1::AbstractVector{GroupIdx1}, coupling_weights1::AbstractMatrix{T} where T<:Real,
    group_ids2::AbstractVector{GroupIdx2}, coupling_weights2::AbstractMatrix{T} where T<:Real,
    JointGroupIdx::DataType=UInt32, PersonIdx::DataType=UInt32) where {GroupIdx1<:Integer, GroupIdx2<:Integer}

  num_groups1 = couplingsizecheck(group_ids1, coupling_weights1)
  num_groups2 = couplingsizecheck(group_ids2, coupling_weights2)

  # grouping2 is the faster changing index
  joint_ids = jointgrouping.(group_ids2, num_groups2, group_ids1) .|> JointGroupIdx
  grouping = PopulationGrouping(joint_ids, num_groups1 * num_groups2, PersonIdx)

  group_sizes = reshape(groupsizes(grouping), num_groups2, num_groups1)
  group1_sizes = vec(sum(group_sizes, dims=1))
  matrix_sampler1 = MatrixAliasSampler(coupling_weights1  .* group1_sizes, GroupIdx1)

  Prod2CouplingSampler(group_ids1, matrix_sampler1, group_ids2, coupling_weights2, grouping)
end

function sample(rng::AbstractRNG, coupling::Prod2CouplingSampler, source_group1_id::Integer, source_group2_id::Integer)
  # first sample the group1
  target_group1_id = sample(rng, coupling.matrix_sampler1, source_group1_id)

  # extract the relevant column for sampling from group2
  target_prob_column2 = view(coupling.weights_matrix2, :, source_group2_id)

  num_groups2 = length(target_prob_column2)
  # compute the sum of weights
  weight_sum = 0.0
  @simd for i in 1:num_groups2
      joint_group_id = jointgrouping(i, num_groups2, target_group1_id)
      weight_sum += target_prob_column2[i] * groupsize(coupling.grouping, joint_group_id)
  end

  target_weight = weight_sum * rand(rng)

  # direct sample the target group 2
  accumulated_weight = 0.0
  target_group2_id = num_groups2
  for i in 1:num_groups2
      joint_group_id = jointgrouping(i, num_groups2, target_group1_id)
      accumulated_weight += target_prob_column2[i] * groupsize(coupling.grouping, joint_group_id)
      if accumulated_weight >= target_weight
          target_group2_id = i
          break
      end
  end

  target_joint_group_id = jointgrouping(target_group2_id, num_groups2, target_group1_id)
  target_group = getgroup(coupling.grouping, target_joint_group_id)

  # sample an individual from the selected group
  rand(rng, target_group)
end