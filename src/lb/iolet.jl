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
    const feq = equilibrium(iolet.density, iolet_μ, cᵢ, wᵢ)[1]
    sim.next_populations[invdir, from...] = feq
end
