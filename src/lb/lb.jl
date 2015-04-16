module LB

export Playground, SingleRelaxationTime, FluidKernel, SingleRelaxationTime
export NashZeroOrderPressure, FluidStreaming, HalfWayBounceBack, ConstantVelocityIOlet
export ConstantPopulationIOlet, density, momentum, velocity, ParabolicVelocityIOlet

using ..Simulation
using ..Playground: Feature, NOTHING
using ..Indices: GridCoords, index, Indexing, gridcoords
import ..Playground.initialize
import ..Indices.neighbor_index
import Base.length

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
function local_kernel(kernel::LocalKernel, sim::Simulation, site::Integer)
    @assert site > 0 && site < length(sim.indexing)
    const from = gridcoords(sim.indexing, site)
    @inbounds const quantities = LocalQuantities(
        typeof(sim).parameters, from, sim.populations[:, site], sim.lattice)
    @inbounds sim.populations[:, site] += (
        collision(kernel.collision, sim.populations[:, site], quantities.feq))
    for direction in 1:length(sim.lattice.weights)
        @inbounds const to = neighbor_index(sim, from, direction)
        const link = to == 0 ? NOTHING: sim.playground[to]
        const streamer = get(kernel.streamers, link, NULLSTREAMER)
        streaming(streamer, quantities, sim, site, to, direction)
    end
end


type Homogeneous{T <: FloatingPoint} <: Initializer
    density::T
    momentum::DenseVector{T}
end
type NullInitializer <: Initializer; end
const NULLINITIALIZER = NullInitializer()

function initialize(init::Homogeneous, sim::Simulation, index::Integer)
    sim.populations[:, index] = equilibrium(sim.lattice, init.density, init.momentum)
end
initialize(init::Homogeneous, sim::Simulation, coords::GridCoords) =
    initialize(init, sim, index(sim.indexing, coords))
initialize(::NullInitializer, sim::Simulation, index::Integer) = nothing


# Both initialize and local_kernel can be used within loops over all sites
for (name, dictionary) in [(:initialize, :initializers), (:local_kernel, :kernels)]
    @eval begin
        function $name(sim::Simulation)
            dic = sim.$dictionary
            for (site, feature) in enumerate(sim.playground)
                if haskey(dic, feature)
                    $name(dic[feature], sim, site)
                end
            end
        end
    end
end
end
