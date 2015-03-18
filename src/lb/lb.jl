module LatticeBoltzmann
using ..Simulation
using ..Playground.Feature
using ..Indices.index

# Base type for all kernally stuff
abstract Kernel
# Local kernel for each lattice type
abstract LocalKernel <: Kernel
# Base type for all collision kernels
abstract Collision <: Kernel
# Base type for all collision kernels
abstract Streaming <: Kernel
# Base type for all initilizers
abstract Initializer <: Kernel

type FluidKernel <: LocalKernel
    # Collision kernel type and data
    collision::Collision
    # Describes how streaming takes place
    streamers :: Dict{Feature, Streaming}
end
type NullKernel <: LocalKernel; end

include("lattice.jl")
include("thermodynamics.jl")
include("collision.jl")
include("streaming.jl")
include("iolet.jl")

local_kernel(kernel::NullKernel, args...; kargs...) = nothing
function local_kernel{T, I}(kernel::FluidKernel, sim::Simulation{T, I}, indices::Vector{I})
    const from = tuple(index(sim.indexing, indices)...)
    const quantities = LocalQuantities{T, I}(indices, sim.populations[:, from...], sim.lattice)
    sim.populations[:, from...] = collision(
        kernel.collision, sim.populations[:, from...], quantities.feq)
    for direction in 1:length(sim.lattice.weights)
      const to_v = index(sim.indexing, indices + sim.lattice.celerities[:, direction])
      const to = tuple(to_v...)
      const link = all(to_v .== 0) ? playground.NOTHING: sim.playground[to...]
      const streamer = get(kernel.streamers, link, NULLSTREAMER)
      streaming(streamer, quantities, sim, from, to, direction)
    end
end


type Homogeneous{T <: FloatingPoint} <: Initializer
    density::T
    momentum::Vector{T}
end

function initialize{T, I}(init::Homogeneous{T}, sim::Simulation{T, I}, indices::Vector{I})
    const from = tuple(index(sim.indexing, indices)...)
    sim.populations[:, from...] = equilibrium(sim.lattice, init.density, init.momentum)
end
end
