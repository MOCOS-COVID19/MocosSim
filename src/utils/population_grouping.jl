using SortingAlgorithms

struct PopulationGrouping{PersonIdx <: Integer}
  person_ids::Vector{PersonIdx}
  group_ptrs::Vector{PersonIdx}
end

function show(io::IO, g::PopulationGrouping)
  num_groups = numgroups(g)
  print(io, "PopulationGrouping with ", numgroups, " groups")
  for group_id in 1:num_groups
    group = getgroup(g, group_id)
    println("Group ", group_id, " of size ", length(group))
  end
end

function make_group_ptrs(group_ids::AbstractArray{T} where T<:Integer, num_groups::Integer, PersonIdx::DataType=UInt32)
  N = length(group_ids)
  # sorting forces individuals of the same group to be together on the list
  sorted_ids = Vector{PersonIdx}(undef, N) .= sortperm(group_ids; alg=RadixSort) #alloc x3

  group_ptrs = zeros(PersonIdx, num_groups + 1) #alloc
  for group_id in group_ids
    group_ptrs[group_id+1] += 1 # size of a group in each entry
  end

  # group i starts at index group_ptrs[i], ends at index group_ptrs[i+1]-1
  group_ptrs[1] = 1
  for i in 1:num_groups
    group_ptrs[i+1] = group_ptrs[i] + group_ptrs[i+1]
  end
  sorted_ids, group_ptrs
end

PopulationGrouping(group_ids, num_groups::Integer, PersonIdx::DataType=UInt32) =
  PopulationGrouping(make_group_ptrs(group_ids, num_groups, PersonIdx)...)

numgroups(g::PopulationGrouping) = length(g.group_ptrs) - 1
getgroup(g::PopulationGrouping, group_id::Integer) = @view g.person_ids[g.group_ptrs[group_id]:(g.group_ptrs[group_id+1]-1)]
