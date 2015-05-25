using LatBo: SandBox, Simulation
using LatBo.Playground: Feature
using LatBo.LB: Streaming, Collision, FluidKernel, LocalKernel
# the following functions are extended for mock types
import LatBo.LB.velocity
import LatBo.LB.collision
import LatBo.LB.streaming
import LatBo.LB.local_kernel


type MockCollision <: Collision
    calls
    args
    MockCollision() = new(0, nothing)
end
function collision(collider::MockCollision, args...)
    collider.calls += 1
    collider.args = deepcopy(args)
    1.5
end

type MockStreaming <: Streaming
    args :: Vector{Any}
    MockStreaming() = new(Any[])
end
function streaming(
        streamer::MockStreaming, quantities::LocalQuantities, sim::Simulation,
        from::Integer, to::Integer, direction::Integer)
    push!(streamer.args, tuple(deepcopy(quantities), sim, from, to, direction))
end

type MockKernel <: LocalKernel
    value
end

local_kernel(kernel::MockKernel, sim::Simulation, site::Integer, args...) =
    (sim.populations[:, site] = kernel.value)

facts("Local fluid kernel") do
    sim = SandBox(:D2Q9, (4, 4))

    const start = index(sim.indexing, [3, 3])
    const cᵢ = sim.lattice.celerities
    const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
    const finish = index(sim.indexing, [3, 3] + cᵢ[:, direction])

    kernel = FluidKernel(
        MockCollision(),
        {
            convert(Feature, 9) => MockStreaming(),
            convert(Feature, 42) => MockStreaming(),
            convert(Feature, 0) => MockStreaming()
        }
    )
    sim.playground[:] = 0
    sim.playground[finish]= 42
    sim.playground[start] = 9
    const fᵢ = 1 + rand(Float64, size(sim.populations, 1))
    const feq = equilibrium(sim.lattice, fᵢ)
    sim.populations[:, start] = fᵢ

    context("Check mock calls") do

        local_kernel(kernel, sim, start)

        @fact kernel.collision.calls => 1
        @fact length(kernel.collision.args) => 2
        @fact kernel.collision.args[1] => roughly(fᵢ)
        @fact kernel.collision.args[2] => roughly(feq)

        @fact length(kernel.streamers[42].args) => 1
        @fact length(kernel.streamers[9].args) => 1
        @fact length(kernel.streamers[0].args) => size(sim.populations, 1) - 2

        @fact kernel.streamers[42].args[1][2] => exactly(sim)
        @fact kernel.streamers[42].args[1][1].feq => roughly(feq)
        @fact kernel.streamers[42].args[1][3:end] => (start, finish, direction)

        @fact kernel.streamers[9].args[1][2] => exactly(sim)
        zerodir = findfirst([all(cᵢ[:, i] .== 0) for i in 1:size(cᵢ, 2)])
        @fact kernel.streamers[9].args[1][3:end] => (start, start, zerodir)
    end

    context("Loop over all sites") do
        sim.playground[:] = 1
        sim.playground[5] = 2
        sim.kernels = {1 => MockKernel(0), 2 => MockKernel(1)}

        sim.populations[:] = -1
        local_kernel(sim)

        @fact any(sim.populations == -1) => false
        @fact sim.populations[:, 5] => roughly(ones(sim.populations[:, 1]))
    end
end
