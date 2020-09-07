using DataStructures

import DataStructures.lt

import Base: push!, pop!, isempty, empty!, sizehint!
export push!, pop!, isempty!, empty!, sizehint!

struct Earlier <: Base.Order.Ordering end

lt(c::Earlier, x::Event, y::Event) = time(x) == time(y) ? kind(x) < kind(y) : time(x) < time(y)

#struct EventQueue
#  immediates::CircularDeque{Event} 
#  laters::BinaryHeap{Event, Earlier}
#  EventQueue(immediate_size::Integer=128) = new(CircularDeque{Event}(immediate_size), BinaryHeap{Event, Earlier}())
#end

#push!(queue::EventQueue, event::Event; immediate::Bool=false) = begin immediate ? push!(queue.immediates, event) : push!(queue.laters, event); queue end
#pop!(queue::EventQueue) = isempty(queue.immediates) ? pop!(queue.laters) : popfirst!(queue.immediates)
#isempty(queue::EventQueue) = isempty(queue.immediates) & isempty(queue.laters)

struct EventQueue
  queue::BinaryHeap{Event, Earlier}
  EventQueue() = new( BinaryHeap{Event, Earlier}())
end

empty!(queue::EventQueue) = empty!(queue.queue.valtree)
push!(queue::EventQueue, event::Event; immediate::Bool=false) = push!(queue.queue, event)
pop!(queue::EventQueue) = pop!(queue.queue)
isempty(queue::EventQueue) = isempty(queue.queue)
sizehint!(queue::EventQueue, sz::Integer) = sizehint!(queue.queue, sz)

