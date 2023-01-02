import Base: iterate, length, eltype, sizehint!, push!

const InfectionIdx = UInt32 # numbers infections, UInt32 -> max 2^32-1 infections
const EdgeIdx = UInt32 # number of edges per infection, UInt32 => max 2^32-1 infected per infection

struct EdgeDictKey
  infection::InfectionIdx
  edgeidx::EdgeIdx
end

struct RobinForest
  infections::Vector{Event} # of length num_infections, for each infection one source Event, in chronological order
  outdegrees::Vector{EdgeIdx} # of length num_infections, for each infection the number of infected

  outedgedict::RobinDict{EdgeDictKey, InfectionIdx} # maps (infection, outedge) to the infection

  recent::Vector{InfectionIdx} # vector of size num_individuals, id of most recent infection for every Person

  RobinForest(num_individuals::Integer) = new(
    Vector{Event}(),
    Vector{EdgeIdx}(),
    RobinDict{EdgeDictKey, InfectionIdx}(),

    fill(InfectionIdx(0), num_individuals)
  )
end

function sizehint!(forest::RobinForest, num_infections::Integer)
  sizehint!(forest.inedges, num_infections)
  sizehint!(forest.outedges, num_infections)
  sizehint!(forest.outedgedict, num_infections)
end

function reset!(forest::RobinForest)
    empty!(forest.infections)
    empty!(forest.outdegrees)
    empty!(forest.outedgedict)
    fill!(forest.recent, InfectionIdx(0))
    nothing
end

function push!(forest::RobinForest, infection::Event)
  push!(forest.infections, infection)
  push!(forest.outdegrees, 0)

  infection_id = length(forest.infections)
  @assert length(forest.outdegrees) == infection_id

  subject_id = subject(infection)
  forest.recent[subject_id] = infection_id

  source_id = source(infection)
  if 0 == source_id
    @assert contactkind(infection) == OutsideContact
    return nothing
  end

  source_infection_id = forest.recent[source_id]
  source_outdeg = forest.outdegrees[source_infection_id] += 1
  outedge_key = EdgeDictKey(source_infection_id, source_outdeg)
  @assert !haskey(forest.outedgedict, outedge_key)
  forest.outedgedict[outedge_key] = infection_id
  nothing
end

recentinfectionof(forest::RobinForest, person_id::Integer) = forest.recent[person_id]
recentbackwardinfectionof(forest::RobinForest, person_id::Integer)::Event = begin
  id = recentinfectionof(forest, person_id)
  0 == id ? Event() : forest.infections[id]
end
recentforwardinfectionsof(forest::RobinForest, person_id::Integer)::OutEdgeList = OutEdgeList(forest, person_id)

recentstrainof(forest::RobinForest, person_id::Integer)::StrainKind = strainkind(recentbackwardinfectionof(forest, person_id))

import Base.iterate
import Base.length
import Base.eltype

struct OutEdgeList
  outedgedict::RobinDict{EdgeDictKey,InfectionIdx}
  infections::Vector{Event}
  infection_id::InfectionIdx
  count::EdgeIdx
  OutEdgeList(forest::RobinForest) = new(forest.outedgedict, forest.inedges, 0, 0)
  OutEdgeList(forest::RobinForest, person_id::Integer) = begin
    infection_id = recentinfectionof(forest, person_id)
    count = 0 == infection_id ? 0 : forest.outdegrees[infection_id]
    new(forest.outedgedict, forest.infections, infection_id, count)
  end
end

iterate(oel::OutEdgeList, state=1) = state > oel.count ? nothing :
    (oel.infections[oel.outedgedict[EdgeDictKey(oel.infection_id, state)]], state+1)
eltype(::Type{OutEdgeList}) = Event
length(oel::OutEdgeList) = oel.count

function saveparams(dict, forest::RobinForest, prefix::AbstractString="")
  N = length(forest.inedges)
  infection_times = Vector{Float32}(undef, N)
  infection_sources = Vector{UInt32}(undef, N)
  contact_kinds = Vector{UInt8}(undef, N)
  strains = Vector{UInt8}(undef, N)
  @simd for i in 1:N
    event = forest.inedges[i]
    infected = NoContact != contactkind(event)
    infection_times[i] = infected ? Float32(time(event)) : NaN32
    infection_sources[i] = infected ? source(event) : 0
    contact_kinds[i] = contactkind(event) |> UInt8
    strains[i] = strainof(forest, i)
  end
  dict[prefix*"infection_times"] = infection_times
  dict[prefix*"infection_sources"] = infection_sources
  dict[prefix*"contact_kinds"] = contact_kinds
  dict[prefix*"strains"] = strains
end