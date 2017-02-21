# Defines streaming types and operators
abstract WallStreaming <: Streaming
abstract IOLetStreaming <: Streaming
immutable FluidStreaming <: Streaming; end
immutable NullStreaming <: Streaming; end
immutable HalfWayBounceBack <: Streaming; end
const NULLSTREAMER = NullStreaming()

# Do nothing when streaming to null, by default
streaming(::NullStreaming, args...) = nothing


# Overload to make unit-testing a bit easier
streaming(
    streamer::Streaming, quant::LocalQuantities, sim::Simulation,
    from::GridCoords, to::GridCoords, direction::Integer
) = streaming(streamer, quant, indexing(sim.indexing, from), indexing(sim.indexing, to), direction)
# Normal streaming from fluid to fluid
streaming(
    streamer::FluidStreaming, quantities::LocalQuantities, sim::Simulation,
    from::Integer, to::Integer, direction::Integer
) = streaming(streamer, sim, from, to, direction)
function streaming(::FluidStreaming, sim::Simulation, from::Integer, to::Integer, dir::Integer)
    # Replacing multi-dimensional getitem with hand-crafted version gains ~10% overall speed-up.
    # At least in this one case. Note that this one case is called very very often.
    const N = size(sim.next_populations, 1)
    @inbounds sim.next_populations[N * (to - 1) + dir] = sim.populations[N * (from - 1) + dir]
end

# Half-way bounce back streaming
streaming(
    streamer::HalfWayBounceBack, quantities::LocalQuantities, sim::Simulation,
    from::Integer, to::Integer, direction::Integer
) = streaming(streamer, sim, from, direction)
function streaming(::HalfWayBounceBack, sim::Simulation, from::Integer, direction::Integer)
    const invdir = sim.lattice.inversion[direction]
    sim.next_populations[invdir, from] = sim.populations[direction, from]
end
