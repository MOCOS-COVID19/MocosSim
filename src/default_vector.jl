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

Base.setindex!(self::DefaultVector{T}, value::T, index::Int64) where T = setindex_impl!(self, value, index)

function setindex_impl!(self::DefaultVector{T}, value::T, index::Int64) where T
    while index > length(self.core)
        push!(self.core, self.defaultValue)
        end
        self.core[index] = value
 end

Base.getindex(self::DefaultVector{T}, index::Int64) where T = self.core[index]

Base.convert(::Type{Vector}, self::DefaultVector{T}) where T = copy(self.core)