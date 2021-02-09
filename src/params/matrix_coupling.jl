function make_group_ptrs(group_ids, num_groups, PersonIdx::DataType=UInt32, GroupIdx::DataType=UInt32)
  N = length(group_ids)
  num_groups = 5
  sorted_ids = sortperm!(Vector{PersonIdx}(undef, N), group_ids)#alloc x2
  counts = zeros(GroupIdx, num_groups) #alloc
  for val in group_ids
    counts[val] += 1
  end
  ranges = fill(UnitRange(GroupIdx(1),GroupIdx(0)), num_groups) #alloc
  ptr = 1
  for i in 1:num_groups
    ranges[i] = UnitRange(ptr, ptr + counts[i] - 1)
    ptr+=counts[i]
  end
  sorted_ids, ranges
end