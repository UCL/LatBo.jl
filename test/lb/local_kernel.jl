using LatBo: SandBox
using LatBo.Playground: Feature
using LatBo.LatticeBoltzmann: Streaming, Collision, FluidKernel, local_kernel
# the following functions are extended for mock types
import LatBo.LatticeBoltzmann.velocity
import LatBo.LatticeBoltzmann.collision
import LatBo.LatticeBoltzmann.streaming


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
function streaming(streamer::MockStreaming, quantities, sim, from, to, direction)
    push!(streamer.args, tuple(deepcopy(quantities), sim, from, to, direction))
end

facts("Local fluid kernel") do
    sim = SandBox(:D2Q9, (4, 4))

    const start = [3, 3]
    const cᵢ = sim.lattice.celerities
    const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
    const finish = start + cᵢ[:, direction]

    kernel = FluidKernel(
        MockCollision(),
        {
            convert(Feature, 9) => MockStreaming(),
            convert(Feature, 42) => MockStreaming(),
            convert(Feature, 0) => MockStreaming()
        }
    )
    sim.playground[:] = 0
    sim.playground[finish...] = 42
    sim.playground[start...] = 9
    const fᵢ = 1 + rand(Float64, size(sim.populations, 1))
    const feq = equilibrium(sim.lattice, fᵢ)
    sim.populations[:, start...] = fᵢ

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
        @fact kernel.streamers[42].args[1][3] => tuple(start...)
        @fact kernel.streamers[42].args[1][4] => tuple(finish...)
        @fact kernel.streamers[42].args[1][5] => direction

        @fact kernel.streamers[9].args[1][2] => exactly(sim)
        @fact kernel.streamers[9].args[1][3] => tuple(start...)
        @fact kernel.streamers[9].args[1][4] => tuple(start...)
        zerodir = findfirst([all(cᵢ[:, i] .== 0) for i in 1:size(cᵢ, 2)])
        @fact kernel.streamers[9].args[1][5] => zerodir
    end
end
