export lbgk
using SIUnits.ShortUnits
using .Units: Time, Length, Viscosity, Density, Velocity, Pressure, Weight, mmHg, LBUnits
using .LB: Lattice

# Sets up calculation with standard lbgk
function lbgk(
        lattice::Lattice, dimensions::(Integer...), δt::Time, δx::Length;
        viscosity::Viscosity=1e-3Pa*s, ρ₀::Density=1e3kg/m^3,
        μ₀::typeof([0., 0] * (m/s))=[0., 0] * (m/s), p₀::Pressure=as(80.0mmHg, Pa), Δp::Pressure=0Pa,
        kwargs...)

    const T = typeof(lattice).parameters[1]
    const units = LBUnits(T, δt, δx, ρ₀)

    const μ = T[u for u in as(units, μ₀)]
    const ρ = convert(T, as(units, ρ₀))
    const p = convert(T, as(units, p₀))
    const c_s² = 1./3.
    const c_s = sqrt(c_s²)
    const τ = 1. + as(units, viscosity) / c_s²

    density(p) = 1. + as(units, p) / c_s²

    # Set up main item
    result = SandBox(lattice, dimensions; ρ₀=ρ, μ₀=μ, kwargs...)
    # Set up kernel and streamers with some default values
    streamers = {
        Playground.INLET  => NashZeroOrderPressure{T}(T[0, 1], ρ),
        Playground.OUTLET => NashZeroOrderPressure{T}(T[0, 1], convert(T, density(p₀ - Δp))),
        Playground.FLUID  => FluidStreaming(),
        Playground.SOLID  => HalfWayBounceBack()
    }
    result.kernels = {Playground.FLUID => FluidKernel(SingleRelaxationTime{T}(1.0/τ), streamers)}

    # Set up playground as box with inflow on left, outflow on right
    if length(dimensions) == 2
        Playground.initialize(result.playground) do i, j
            if j == 1 || j == dimensions[2]
                return Playground.SOLID
            elseif i == 1
                return Playground.INLET
            elseif i == dimensions[1]
                return Playground.OUTLET
            else
                return Playground.FLUID
            end
        end
    else
        Playground.initialize(result.playground) do i, j, k
            if j == 1 || j == dimensions[2] || k == 1 || k == dimensions[3]
                return Playground.SOLID
            elseif i == 1
                return Playground.INLET
            elseif i == dimensions[1]
                return Playground.OUTLET
            else
                return Playground.FLUID
            end
        end
    end
    result
end

function lbgk(lattice::Symbol, dimensions::(Integer...), δt::Time, δx::Length; kwargs...)
    lbgk(getfield(LB, lattice)::LB.Lattice, dimensions, δt, δx; kwargs...)
end
