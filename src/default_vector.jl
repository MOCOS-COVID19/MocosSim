#=
deafult_vector:
- Julia version: 
- Author: jacek
- Date: 2022-03-18
=#

mutable struct DefaultVector{T}
    defaultValue::T
    core::Vector{T}

    DefaultVector(defaultValue_) = new{typeof(defaultValue_)}(defaultValue_, Vector{typeof(defaultValue_)}(undef, 0))
    DefaultVector(defaultValue_, K) = new{K}(defaultValue_, Vector{K}(undef, 0))
end

Base.setindex!(self::DefaultVector{T}, value::T, index::Int64) where T = _ensure_index_in_range!(self, index) do
    self.core[index] = value
end

Base.getindex(self::DefaultVector{T}, index::Int64) where T = _ensure_index_in_range!(self, index) do
    self.core[index]
end

Base.convert(::Type{Vector}, self::DefaultVector{T}) where T = copy(self.core)

function _ensure_index_in_range!(next_func::Function, self::DefaultVector{T}, index::Int64) where T
    while index > length(self.core)
        push!(self.core, self.defaultValue)
    end
    next_func()
end