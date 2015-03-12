using FactCheck: facts, @fact, roughly
using LatBo: kernel
using LatBo.kernel: collision, index, Cartesian, Periodic

facts("Kernel actions") do
    context("SRT collision kernel") do
        # random number we can play with and avoid overflow
        τ⁻¹ = convert(Int64, 2)
        f = convert(Array{Int64}, rand(Int8, 3))
        @fact collision(τ⁻¹, f, f) => zeros(Int16, length(f))
        @fact collision(τ⁻¹*τ⁻¹, f, f) => zeros(Int16, length(f))
        @fact collision(τ⁻¹, 2*f, f) => -τ⁻¹ * f
        @fact collision(τ⁻¹, f, 2*f) => τ⁻¹ * f
        @fact collision(τ⁻¹, 4*f, 2*f) => collision(2*τ⁻¹, 2*f, f)
    end

    context("Cartesian indexing") do
        @fact index(Cartesian(), [1 2 3]) => [1 2 3]
        @fact index(Cartesian(), [8 10 20]) => [8 10 20]
        @fact index(Cartesian(), "hello") => "hello"
    end

    context("Periodic indexing") do
        @fact index(Periodic((64, 64, 64)), [8, 16, 32]) => [8, 16, 32]
        @fact index(Periodic((5, 64, 64)), [8, 16, 32]) => [3, 16, 32]
        # Between 1 and n, not 0 and n-1. Julia indices start at 1
        @fact index(Periodic((64, 16, 64)), [8, 0, 32]) => [8, 16, 32]
        @fact index(Periodic((64, 64, 32)), [8, 16, 32]) => [8, 16, 32]
        @fact index(Periodic((64, 64, 32)), [8, 16, 33]) => [8, 16, 1]
    end
end
