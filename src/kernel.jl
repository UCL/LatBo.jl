module lb
using ..Collision
using ..Streaming
using ..Indexing
using ..thermodynamics
using ..speed_of_sound_squared
using ..Feature
using ..LocalKernel
using ..Simulation
using ..thermodynamics.LocalQuantities

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


# Normal streaming from fluid to fluid
streaming{T, I}(
    streamer::FluidStreaming, quantities::LocalQuantities{T, I}, sim::Simulation{T, I},
    from::(I...), to::(I...), direction::I
) = streaming(streamer, from, to, direction)
function streaming{T, I}(
    ::FluidStreaming, sim::Simulation{T, I}, from::(I...), to::(I...), direction::I)

    @assert(length(from) == length(to))
    sim.next_populations[direction, to...] = sim.populations[direction, from...]
end

# Half-way bounce back streaming
streaming{T, I}(
    streamer::HalfWayBounceBack, quantities::LocalQuantities{T, I}, sim::Simulation{T, I},
    from::(I...), to::(I...), direction::I
) = streaming(streamer, sim, from, direction)
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
    maxspeed :: T

    function ParabolicVelocityIOlet(n₀, Γ, r₀, ν_max)
        n = T[n₀...]
        @assert ndims(n₀) == 1 && ndims(Γ) == 1 && length(Γ) == length(n₀)
        @assert dot(n, n) > 1e-12 && r₀ > 1e-12 && ν_max > 1e-12
        new(n / norm(n), T[Γ...], r₀, ν_max)
    end
end

streaming{T, I}(
    streamer::VelocityIOlet, quantities::LocalQuantities{T, I}, sim::Simulation{T, I},
    from::(I...), to::(I...), direction::I
) = streaming(streamer, sim, from, to, direction)
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

# Time-independent parabolic inlet
velocity(iolet::ParabolicVelocityIOlet, position, time) = velocity(iolet, position)
function velocity(iolet::ParabolicVelocityIOlet, position)
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
    function NashZeroOrderPressure(n₀, ρ)
        n = T[n₀...]
        @assert dot(n, n) > 1e-12 && ρ > 1e-12
        new(n / norm(n), ρ)
    end
end

streaming{T, I}(
    iolet::NashZeroOrderPressure, quantitie::LocalQuantities{T, I}, sim::Simulation{T, I},
    from::(I...), to::(I...), direction::I
) = streaming(iolet, quantities.velocity, sim, from, direction)
function streaming{T, I}(
    iolet::NashZeroOrderPressure, velocity::Vector{T}, sim::Simulation{T, I},
    from::(I...), direction::I)
    const invdir = sim.lattice.inversion[direction]
    const iolet_μ = iolet.density * dot(iolet.normal, velocity) * iolet.normal
    const cᵢ = sim.lattice.celerities[:, invdir:invdir]
    const wᵢ = sim.lattice.weights[invdir:invdir]
    const feq = thermodynamics.equilibrium(iolet.density, iolet_μ, cᵢ, wᵢ)[1]
    sim.next_populations[invdir, from...] = feq
end

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

end
