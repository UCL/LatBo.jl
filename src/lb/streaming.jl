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
