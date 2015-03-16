module lb
using LatBo: Collision, Streaming, Indexing, thermodynamics, speed_of_sound_squared, Feature,
        LocalKernel, Simulation
using LatBo.thermodynamics: LocalQuantities

type FluidKernel <: LocalKernel
    # Collision kernel type and data
    collision::Collision
    # Describes how streaming takes place
    streamers :: Dict{Feature, Streaming}
end
type NullKernel <: LocalKernel; end

type SingleRelaxationTime{T <: Real} <: Collision
    # Inverse of the relaxation time
    τ⁻¹::T
end

# Defines collision kernels for SRT
collision{T}(τ⁻¹::T, fᵢ::Vector{T}, feq::Vector{T}) = τ⁻¹ * (feq - fᵢ)
collision{T}(k::SingleRelaxationTime{T}, fᵢ::Vector{T}, feq::Vector{T}) = collision(k.τ⁻¹, feq, fᵢ)

# No frills indexing
immutable Cartesian <: Indexing
    dimensions::Vector{Int64}
end
# Indices are periodic
immutable Periodic <: Indexing
    dimensions::Vector{Int64}
end

# A function to retrieve array indices from simulation indices
# Dumps to zero if indices are out of bounds
function index{T1 <: Int, T2 <: Int}(dimensions::Vector{T1}, indices::Vector{T2})
    @assert(size(dimensions) == size(indices))
    any(indices .< 1) || any(indices .> dimensions) ?
        zeros(T2, size(indices)): indices
end
index(kernel::Indexing, indices) = indices
index(kernel::Cartesian, indices) = index(kernel.dimensions, indices)
index{T <: Int}(kernel::Periodic, indices::Array{T}) = 1 + mod(indices - 1, kernel.dimensions)

# Defines streaming types and operators
abstract WallStreaming <: Streaming
abstract IOLetStreaming <: Streaming
immutable type FluidStreaming <: Streaming; end
immutable type NullStreaming <: Streaming; end
immutable type HalfWayBounceBack <: Streaming; end
const NULLSTREAMER = NullStreaming()

# Do nothing when streaming to null, by default
streaming(::NullStreaming, args...) = nothing
# Some overloading to simplify specialized functions where possible
streaming{T, I}(
    streamer::Streaming, quantities::LocalQuantities{T, I}, sim::Simulation{T, I},
    from::(I...), to::(I...), direction::I
) = streaming(streamer, sim, from, to, direction)
streaming{T, I}(
    streamer::Streaming, quantities::LocalQuantities{T, I}, sim::Simulation{T, I},
    from::(I...), to::(I...), direction::I
) = streaming(streamer, quantities, sim, from, direction)
streaming{T, I}(
    streamer::Streaming, sim::Simulation{T, I}, from::(I...), to::(I...), direction::I
) = streaming{T, I}(streamer, sim, from, direction)
# Normal streaming from fluid to fluid
function streaming{T, I}(
    ::FluidStreaming, sim::Simulation{T, I}, from::(I...), to::(I...), direction::I)

    @assert(length(from) == length(to))
    sim.next_populations[direction, to...] = sim.populations[direction, from...]
end
# Half-way bounce back streaming
function streaming{T, I}(::HalfWayBounceBack, sim::Simulation{T, I}, from::(I...), direction::I)
    const invdir = sim.lattice.inversion[direction]
    sim.next_populations[invdir, from...] = sim.populations[direction, from...]
end


# Create a parabolic velocity inlet
abstract VelocityIOlet <: IOLetStreaming
type ParabolicVelocityIOlet{T} <: VelocityIOlet
    # Direction of the inlet
    normal :: Vector{T}
    # Center of the inlet
    center :: Vector{T}
    # Radius of the inlet
    radius :: T
    # Maximum speed
    maxspeed :: Vector{T}
end

function streaming{T, I}(
        streamer::VelocityIOlet, sim::Simulation{T, I}, from::(I...), to::(I...), direction::I)
    const bb_direction = sim.lattice.inversion[direction]
    const cᵢ = sim.lattice.celerities[:, direction]
    const wᵢ = sim.lattice.weights[direction]
    const halfway = T[from...] + 0.5cᵢ
    const μ = velocity(streamer, halfway, sim.time)
    const correction = 2wᵢ * dot(cᵢ, μ) / speed_of_sound_squared
    sim.next_populations[bb_direction, from...] = sim.populations[direction, from...] - correction
end

function velocity(iolet::ParabolicVelocityIOlet, position, time)
    const centered = position - iolet.center
    const z = dot(centered, iolet.normal)
    const radial = (dot(centered, centered) - z * z) / iolet.radius / iolet.radius
    ((1 - radial) * iolet.maxspeed) * iolet.normal
end


abstract PressureIOlet{T} <: IOLetStreaming
type NashZeroOrderPressure{T} <: IOLetStreaming
    # Direction of the inlet
    normal :: Vector{T}
    # Density at the inlet
    density :: T
end

function streaming{T, I}(
    iolet::NashZeroOrderPressure, quantitie::LocalQuantities{T, I}, sim::Simulation{T, I},
    from::(I...), direction::I)
    const invdir = sim.lattice.inversion[direction]
    const iolet_μ = iolet.density * dot(iolet.normal, quantities.velocity) * iolet.normal
    const cᵢ = sim.lattice.celerities[:, invdir:invdir]
    const wᵢ = sim.lattice.weights[invdir:invdir]
    const feq = thermodynamics.equilibrium(iolet.density, iolet_μ, cᵢ, wᵢ)[1]
    sim.next_populations[invdir, from...] = feq
end

local_kernel(kernel::NullKernel, args...; kargs...) = nothing
function local_kernel{T, I}(kernel::FluidKernel, sim::Simulation{T, I}, indices::Vector{I})
    const from = indexing(sim.indexing, indices)
    const quantities = LocalQuantities(indices, sim.population[:, from], sim.lattice)
    sim.populations[:, from] = collision(kernel.collision, sim.population[:, from], quantities.feq)
    for direction in 1:length(sim.lattice.weights)
      const to_v = index(sim.indexing, indices + sim.lattice.celerities[direction])
      const to = to_v...
      const link = all(to_v .== 0) ? playground.NOTHING: sim.playground[to]
      const streamer = get(kernel.streamers, link, NULLSTREAMER)
      streaming(streamer, quantities, sim, from, to, direction)
    end
end

end
