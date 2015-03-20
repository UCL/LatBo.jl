module LB
using ..Simulation
using ..Playground.Feature
using ..Indices.GridCoords
using ..Indices.index
using ..Indices.gridcoords

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
function local_kernel(kernel::FluidKernel, sim::Simulation, site::Integer)
    const from = gridcoords(sim.indexing, site)
    const quantities = LocalQuantities(
        typeof(sim).parameters, from, sim.populations[:, site], sim.lattice)
    sim.populations[:, site] = collision(kernel.collision, sim.populations[:, site], quantities.feq)
    for direction in 1:length(sim.lattice.weights)
        const to = index(sim.indexing, from + sim.lattice.celerities[:, direction])
        const link = to == 0 ? playground.NOTHING: sim.playground[to]
        const streamer = get(kernel.streamers, link, NULLSTREAMER)
        streaming(streamer, quantities, sim, site, to, direction)
    end
end


type Homogeneous{T <: FloatingPoint} <: Initializer
    density::T
    momentum::Vector{T}
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
# So run them both together
for (name, dictionary) in [(:initialize, :initializers), (:local_kernel, :kernels)]
    @eval begin
        function $name(sim::Simulation)
            dic = sim.$dictionary
            for (site, feature) in enumerate(sim.playground)
                $name(get(dic, feature, NULLINITIALIZER), sim, site)
            end
        end
    end
end

end
