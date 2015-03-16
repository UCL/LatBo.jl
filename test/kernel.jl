using FactCheck: facts, @fact, roughly
using LatBo: Indexing, D3Q19, Simulation, speed_of_sound_squared
using LatBo.lb: collision, index, Cartesian, Periodic, streaming, FluidStreaming,
        HalfWayBounceBack, IOLetStreaming, VelocityIOlet
import LatBo.lb.velocity

type DummyVelocityIOlet <: VelocityIOlet
    position::Vector{Float64}
    time::Float64
    DummyVelocityIOlet() = new(Float64[], 0)
end
function velocity(streamer::DummyVelocityIOlet, position, time)
    streamer.position = deepcopy(position)
    streamer.time = time
    Float64[-10, 10]
end


facts("Kernel actions") do
    context("Cartesian indexing") do
        @fact index(Cartesian([10, 10, 10]), [1, 2, 3]) => [1, 2, 3]
        @fact index(Cartesian([10, 10, 30]), [8, 10, 20]) => [8, 10, 20]
        @fact index(Cartesian([10, 10, 30]), [0, 10, 20]) => [0, 0, 0]
        @fact index(Cartesian([10, 10, 30]), [10, 11, 20]) => [0, 0, 0]
    end

    context("Periodic indexing") do
        @fact index(Periodic([64, 64, 64]), [8, 16, 32]) => [8, 16, 32]
        @fact index(Periodic([5, 64, 64]), [8, 16, 32]) => [3, 16, 32]
        # Between 1 and n, not 0 and n-1. Julia indices start at 1
        @fact index(Periodic([64, 16, 64]), [8, 0, 32]) => [8, 16, 32]
        @fact index(Periodic([64, 64, 32]), [8, 16, 32]) => [8, 16, 32]
        @fact index(Periodic([64, 64, 32]), [8, 16, 33]) => [8, 16, 1]
    end

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

    sim = Simulation(D2Q9, (4, 4))

    context("Fluid to Fluid streaming") do
        sim.populations[:] = 0
        sim.next_populations[:] = 0
        const cᵢ = sim.lattice.celerities
        const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
        const start = (3, 3)
        const finish = (2, 4)
        sim.populations[direction, start...] = 1
        @fact sim.next_populations .== 0 => all
        # Perform streaming
        streaming(FluidStreaming(), sim, start, finish, direction)
        # Check it went to the right place and only there
        @fact sim.next_populations[direction, finish...] => roughly(1)
        @fact sum(abs(sim.next_populations)) => roughly(1)
        # Check populations was not affected
        @fact sim.populations[direction, start...] => roughly(1)
        @fact sum(abs(sim.populations)) => roughly(1)
    end

    context("Halfway bounce-back streaming") do
        sim.populations[:] = 0
        sim.next_populations[:] = 0
        const cᵢ = sim.lattice.celerities
        const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
        const invdir = find([all(cᵢ[:, i] .== [1, -1]) for i in 1:size(cᵢ, 2)])[1]
        const start = (3, 3)

        @fact sim.lattice.inversion[direction] => invdir
        sim.populations[direction, start...] = 1
        @fact sim.next_populations .== 0 => all
        # Perform streaming
        streaming(HalfWayBounceBack(), sim, start, direction)
        # Check it went to the right place and only there
        @fact sim.next_populations[invdir, start...] => roughly(1)
        @fact sum(abs(sim.next_populations)) => roughly(1)
        @fact sim.populations[direction, start...] => roughly(1)
        @fact sum(abs(sim.populations)) => roughly(1)
    end

    context("Velocity iolet") do

        const cᵢ = sim.lattice.celerities
        const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
        const invdir = find([all(cᵢ[:, i] .== [1, -1]) for i in 1:size(cᵢ, 2)])[1]
        const start = (3, 3)
        const finish = (2, 4)
        const halfway = Float64[start...] + 0.5*Float64[-1, 1]
        sim.time = 40.0


        sim.populations[:] = 0
        sim.next_populations[:] = 0
        @fact sim.lattice.inversion[direction] => invdir
        sim.populations[direction, start...] = 1
        @fact sim.next_populations .== 0 => all

        # perform streaming
        iolet = DummyVelocityIOlet()
        streaming(iolet, sim, start, finish, direction)

        @fact iolet.position => roughly(Float64[2.5, 3.5])
        @fact iolet.time => roughly(sim.time)
        # 20 is the celerity * momentum
        expected = 2sim.lattice.weights[direction] / speed_of_sound_squared * 20
        @fact sim.next_populations[invdir, start...] => roughly(1 - expected)
        @fact sum(abs(sim.next_populations)) => roughly(abs(1 - expected))
        @fact sim.populations[direction, start...] => roughly(1)
        @fact sum(abs(sim.populations)) => roughly(1)
    end
end
