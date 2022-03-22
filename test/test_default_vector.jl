#=
test_default_vector:
- Julia version: 
- Author: jacek
- Date: 2022-03-18
=#
using Test

include("../src/default_vector.jl")


function test_default_vector()
    for v in [DefaultVector(1), DefaultVector(1, Int)]
        v[1] = 5
        v[3] = 7
        @test v[2] == 1
        @test convert(Vector, v) == [5, 1, 7]

        v[7] += 1

        @test convert(Vector, v) == [5, 1, 7, 1, 1, 1, 2]
    end

    for v in [DefaultVector('c'), DefaultVector('c', typeof('c'))]
        v[3] = 'a'
        @test convert(Vector, v) == ['c', 'c', 'a']
    end
end

test_default_vector()