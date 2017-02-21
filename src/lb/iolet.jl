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

streaming(
    streamer::VelocityIOlet, quantities::LocalQuantities, sim::Simulation,
    from::Integer, to::Integer, direction::Integer
) = streaming(streamer, sim, from, direction)
function streaming(streamer::VelocityIOlet, sim::Simulation, from::Integer, direction::Integer)
    const bb_direction = sim.lattice.inversion[direction]
    const cᵢ = sim.lattice.celerities[:, direction]
    const wᵢ = sim.lattice.weights[direction]
    const halfway = gridcoords(sim.indexing, from) + 0.5cᵢ
    const μ = velocity(streamer, halfway, sim.time)
    const correction = 2wᵢ * dot(cᵢ, μ) / speed_of_sound_squared
    sim.next_populations[bb_direction, from] = sim.populations[direction, from] - correction
end

# Time-independent parabolic inlet
velocity(iolet::ParabolicVelocityIOlet, position, time) = velocity(iolet, position)
function velocity(iolet::ParabolicVelocityIOlet, position)
    const centered = position - iolet.center
    const z = dot(centered, iolet.normal)
    const radial = (dot(centered, centered) - z * z) / iolet.radius / iolet.radius
    ((1 - radial) * iolet.maxspeed) * iolet.normal
end

type ConstantVelocityIOlet{T <: AbstractFloat} <: VelocityIOlet
    velocity :: Vector{T}
end
streaming(
    streamer::ConstantVelocityIOlet, quantities::LocalQuantities, sim::Simulation,
    from::Integer, to::Integer, direction::Integer
) = streaming(streamer, sim, from, direction)
function streaming(iolet::ConstantVelocityIOlet, sim::Simulation, from::Integer, direction::Integer)
    const bb_direction = sim.lattice.inversion[direction]
    const cᵢ = sim.lattice.celerities[:, direction]
    const wᵢ = sim.lattice.weights[direction]
    const correction = 2wᵢ * dot(cᵢ, iolet.velocity) / speed_of_sound_squared
    sim.next_populations[bb_direction, from] = sim.populations[direction, from] - correction
end


type NashZeroOrderPressure{T} <: IOLetStreaming
    # Direction of the inlet
    normal :: Vector{T}
    # Density at the inlet
    density :: T
    function NashZeroOrderPressure(n₀::Vector, ρ::Number)
        n = T[n₀...]
        @assert dot(n, n) > 1e-12 && ρ > 1e-12
        new(n / norm(n), ρ)
    end
end

streaming(
    iolet::NashZeroOrderPressure, quantities::LocalQuantities, sim::Simulation,
    from::Integer, to::Integer, direction::Integer
) = streaming(iolet, quantities.velocity, sim, from, direction)
function streaming(
        iolet::NashZeroOrderPressure, velocity::Vector, sim::Simulation,
        from::Integer, direction::Integer)
    const invdir = sim.lattice.inversion[direction]
    const iolet_μ = iolet.density * dot(iolet.normal, velocity) * iolet.normal
    const cᵢ = sim.lattice.celerities[:, invdir:invdir]
    const wᵢ = sim.lattice.weights[invdir:invdir]
    const feq = equilibrium(iolet.density, iolet_μ, cᵢ, wᵢ)[1]
    sim.next_populations[invdir, from] = feq
end

type ConstantPopulationIOlet{I <: Integer} <: IOLetStreaming
    # Normal to iolet
    # Should be a lattice vector
    # Should points from fluid site to iolet site
    normal :: Vector{I}
end

function streaming(
    iolet::ConstantPopulationIOlet, quantities::LocalQuantities, sim::Simulation,
    from::Integer, to::Integer, direction::Integer)
    # This object can stream to a different object than expected
    const invdir = sim.lattice.inversion[direction]
    const mirror = gridcoords(sim.indexing, from) + iolet.normal + sim.lattice.celerities[:, invdir]
    const mirror_site = index(sim.indexing, mirror)
    const link = to == 0 ? NOTHING: sim.playground[mirror_site]
    const this_kernel = sim.kernels[sim.playground[from]]
    const streamer = get(this_kernel.streamers, link, NULLSTREAMER)
    @assert !(typeof(streamer) <: IOLetStreaming)
    streaming(streamer, quantities, sim, from, mirror_site, invdir)
end
