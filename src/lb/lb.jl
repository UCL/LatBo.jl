module LB

export Playground, FluidKernel, SingleRelaxationTime
export NashZeroOrderPressure, FluidStreaming, HalfWayBounceBack, ConstantVelocityIOlet
export ConstantPopulationIOlet, density, momentum, velocity, ParabolicVelocityIOlet

using ..AbstractLattice
using ..Simulation
using ..Playground: Feature, NOTHING
using ..Indices: GridCoords, index, Indexing, gridcoords
import ..Playground.initialize
import ..Indices.neighbor_index
import Base: length, ndims

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
# Holds local quantities
# - density
# - momentum
# - velocity
# - feq
abstract AbstractLocalQuantities

type FluidKernel <: LocalKernel
    # Collision kernel type and data
    collision::Collision
    # Describes how streaming takes place
    streamers :: Dict{Feature, Streaming}
end
type NullKernel <: LocalKernel; end

include("lattice.jl")
include("thermodynamics.jl")
include("LocalQuantities.jl")
include("collision.jl")
include("streaming.jl")
include("iolet.jl")

local_kernel(kernel::NullKernel, args...; kargs...) = nothing
function local_kernel(kernel::LocalKernel, sim::Simulation, site::Integer)
    quantities = LocalQuantities(sim.lattice)
    fᵢ = similar(sim.lattice.weights)
    local_kernel(kernel, sim, site, quantities, fᵢ)
end

function local_kernel(
    kernel::LocalKernel, sim::Simulation, site::Integer, quantities::LocalQuantities,
    fᵢ::DenseVector)
    @assert site > 0 && site < length(sim.indexing)
    @inbounds fᵢ[:] = sim.populations[:, site]
    @inbounds LocalQuantities!(quantities, fᵢ, sim.lattice)
    @inbounds collision!(sim.populations, kernel.collision, quantities.feq, site)
    for direction in 1:length(sim.lattice.weights)
        @inbounds const to = neighbor_index(sim, site, direction)
        const link = to == 0 ? NOTHING: sim.playground[to]
        const streamer = get(kernel.streamers, link, NULLSTREAMER)
        streaming(streamer, quantities, sim, site, to, direction)
    end
end


type Homogeneous{T <: AbstractFloat} <: Initializer
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


# Loop over all fluid sites and initialize populations
function initialize(sim::Simulation)
    for (site, feature) in enumerate(sim.playground)
        if haskey(sim.initializers, feature)
            initialize(sim.initializers[feature], sim, site)
        end
    end
end
# Loop over all fluid sites and perform calculations
function local_kernel(sim::Simulation)
    quants = LocalQuantities(sim.lattice)
    fᵢ = similar(sim.lattice.weights)
    for (site, feature) in enumerate(sim.playground)
        if haskey(sim.kernels, feature)
            local_kernel(sim.kernels[feature], sim, site, quants, fᵢ)
        end
    end
end
end
