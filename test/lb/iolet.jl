using LatBo: SandBox
using LatBo.LB: VelocityIOlet, ParabolicVelocityIOlet,
        equilibrium, NashZeroOrderPressure, streaming, speed_of_sound_squared
import LatBo.LB: velocity

type MockVelocityIOlet <: VelocityIOlet
    position::Vector{Float64}
    time::Float64
    MockVelocityIOlet() = new(Float64[], 0)
end
function velocity(streamer::MockVelocityIOlet, position, time)
    streamer.position = deepcopy(position)
    streamer.time = time
    Float64[-10, 10]
end

facts("Inlets and Outlets") do
    sim = SandBox(:D2Q9, (4, 4))

    context("Velocity") do
        context("General algorithm") do
            const cᵢ = sim.lattice.celerities
            const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
            const invdir = find([all(cᵢ[:, i] .== [1, -1]) for i in 1:size(cᵢ, 2)])[1]
            const start = index(sim.indexing, (3, 3))
            const finish = index(sim.indexing, (2, 4))
            const halfway = Float64[3, 3] + 0.5*Float64[-1, 1]
            sim.time = 40.0

            sim.populations[:] = 0
            sim.next_populations[:] = 0
            @fact sim.lattice.inversion[direction] --> invdir
            sim.populations[direction, start] = 1
            @fact sim.next_populations .== 0 --> all

            # perform streaming
            iolet = MockVelocityIOlet()
            streaming(iolet, sim, start, direction)

            @fact iolet.position --> roughly(Float64[2.5, 3.5])
            @fact iolet.time --> roughly(sim.time)
            # 20 is the celerity * momentum
            expected = 2sim.lattice.weights[direction] / speed_of_sound_squared * 20
            @fact sim.next_populations[invdir, start] --> roughly(1 - expected)
            @fact sum(abs(sim.next_populations)) --> roughly(abs(1 - expected))
            @fact sim.populations[direction, start] --> roughly(1)
            @fact sum(abs(sim.populations)) --> roughly(1)
        end

        context("parabolic") do
            const n₀, n₁, Γ, r, ν_max = [1., 0], [0., 1], [5., 5], 5., 10.
            const ϵ = 1e-6
            iolet = ParabolicVelocityIOlet{Float64}(n₀, Γ, r, ν_max)

            # should be zero at radius
            @fact velocity(iolet, Γ + r * n₁ + rand() * n₀) --> roughly(zeros(Float64, 2))
            # should be ν_max at origin
            @fact velocity(iolet, Γ) --> roughly(ν_max * n₀)
            # should be 0.75 * ν_max at origin + 1/2 radius
            @fact velocity(iolet, Γ + 0.5r*n₁) --> roughly(0.75ν_max * n₀)

            # velocity should be along normal
            @fact dot(velocity(iolet, Γ + 2r*rand(Float64, 2)), n₁) --> roughly(0)
            # velocity should not depend on position along normal
            α = Γ + 2r*rand(Float64, 2)
            @fact velocity(iolet, α + r * rand() * n₀) --> roughly(velocity(iolet, α))

            # velocity should be maximum at origin and equal to maxspeed
            α = Γ + 2r*rand() * n₀
            δν = velocity(iolet, α + ϵ * n₁) - velocity(iolet, α - ϵ * n₁)
            @fact dot(δν, n₀) --> roughly(0)
        end
    end

    context("Nash zero pressure") do
        const n₀, n₁, ρ = [1., 0], [0., 1], 5.
        const cᵢ = sim.lattice.celerities
        const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
        const invdir = find([all(cᵢ[:, i] .== [1, -1]) for i in 1:size(cᵢ, 2)])[1]
        const start = index(sim.indexing, (3, 3))
        const ν₀ = 5n₀ + 6n₁

        iolet = NashZeroOrderPressure{Float64}(n₀, ρ)

        # expected values
        μ = ρ * dot(ν₀, n₀)n₀
        feq = equilibrium(ρ, μ, sim.lattice.celerities, sim.lattice.weights)

        streaming(iolet, ν₀, sim, start, direction)
        @fact sim.next_populations[invdir, start] --> roughly(feq[invdir])
        @fact sum(abs(sim.next_populations[invdir, start])) --> roughly(feq[invdir])

        # original population has no effect
        sim.populations = 1.0 + rand(Float64, size(sim.populations))
        streaming(iolet, ν₀, sim, start, direction)
        @fact sim.next_populations[invdir, start] --> roughly(feq[invdir])
        @fact sum(abs(sim.next_populations[invdir, start])) --> roughly(feq[invdir])
    end
end

