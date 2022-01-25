
struct SplitAgeCouplingParams
  coupling::Prod2CouplingSampler{UInt32, UInt8, PersonIdx}
end

function SplitAgeCouplingParams(
  ages::AbstractVector{T} where T<:Real,
  genders::Union{Nothing, AbstractVector{Bool}},
  age_thresholds::AbstractVector{T} where T<:Real,
  age_coupling_weights::AbstractMatrix{Float64},
  split_group_ids::AbstractVector{T} where T<:Integer,
  split_coupling_weights::AbstractMatrix{Float64})

  @assert age_thresholds[1] == 0
  num_groups = length(age_thresholds)
  @assert (num_groups,num_groups) == size(coupling_weights)

  age_group_ids = agegroup.((age_thresholds,), ages)
  @assert minimum(age_group_ids) > 0
  if genders !== nothing
    @assert length(ages) == length(genders)
    age_group_ids .= age_group_ids .*2 .+ genders .- 1
  end
  @assert maximum(age_group_ids) <= num_groups <= typemax(GroupIdx) "minmax: $(extrema(age_group_ids)), num=$num_groups, typemax=$(typemax(GroupIdx))"

  coupling_sampler = Prod2CouplingSampler(
    split_group_ids, split_coupling_weights,
    age_group_ids, age_coupling_weights,
    UInt32, PersonIdx)

  SplitAgeCouplingParams(coupling_sampler)
end
