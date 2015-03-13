using LatBo: playground
module kernel

# Base type for all kernally stuff
abstract Kernel
# Base type for all collision kernels
abstract Collision <: Kernel
type SingleRelaxationTime{T <: Real} <: Collision
    # Inverse of the relaxation time
    τ⁻¹::T
end

# Defines collision kernels for SRT
collision{T}(τ⁻¹::T, fᵢ::Vector{T}, feq::Vector{T}) = τ⁻¹ * (feq - fᵢ)
collision{T}(k::SingleRelaxationTime{T}, fᵢ::Vector{T}, feq::Vector{T}) = collision(k.τ⁻¹, feq, fᵢ)

# Base type for all indexing kernels
abstract Indexing <: Kernel
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
index(kernel::Indexing, indices) = index(kernel.dimensions, indices)
index{T <: Int}(kernel::Periodic, indices::Array{T}) = 1 + mod(indices - 1, kernel.dimensions)

abstract Streaming <: Kernel
abstract WallStreaming <: Streaming
abstract IOLetStreaming <: Streaming
type BulkStreaming <: Streaming
end
type NullStreaming <: Streaming
end

function streaming{T <: Int}(simulation, indices::Vector{T}, direction::T)
  const from = index(simulation.indexing, indices)...
  const to_v = index(simulation.indexing, indices + simulation.lattice.celerities[direction])
  const to = to_v...

  const link = all(to_v .== 0) ? playground.NOTHING: simulation.playground[to]
  stream_to_nothing(simulation.streamers[link], simulation, from, to, direction)
end

# Do nothing when streaming to null, by default
streaming(streamer::NullStreaming, args...) = nothing
# Normal bulk streaming
function streaming{T <: Int}(
    streamer::BulkStreaming, simulation, from::(T...), to::(T...), direction::T)

    @assert(length(from) == length(to))
    simulation.next_populations[tuple(direction, to...)] =
        simulation.next_populations[tuple(direction, from...)]
end
# Half-way bounce back streaming
end
