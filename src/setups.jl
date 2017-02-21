export lbgk
using .LB: Lattice
using .Units: dimensionless, LBUnits

# Sets up calculation with standard lbgk
# All units should be in SI: 
#  - time -> s
#  - length -> m
#  - pressure -> Pa
#  - viscosity -> Pa * s
function lbgk{N, I <: Integer}(
        lattice::Lattice, dimensions::NTuple{N, I},
        δt::Number,                       # s
        δx::Number;                       # m
        viscosity::Number=1e-3,           # Pa*s
        μ₀::Vector=[0., 0],               # m/s,
        p₀::Number=80.0 * 133.3223684211, # mmHG * x = Pa
        Δp::Number=0,                     # Pa
        kwargs...)

    const T = typeof(lattice).parameters[1]
    const ρ₀=1e3                    # kg*m^-3
    const δm = ρ₀ * δx^3
    const units = LBUnits{T}(δt, δx, δm)

    const μ = dimensionless(units, :velocity, μ₀)
    const p = dimensionless(units, :pressure, p₀)
    const c_s² = 1./3.
    const c_s = sqrt(c_s²)
    const τ = 1. + dimensionless(units, :viscosity, viscosity) / c_s²

    density(p) = 1. + dimensionless(units, :pressure, p) / c_s²

    # Set up main item
    result = SandBox(lattice, dimensions; ρ₀=density(p₀ - 0.5Δp), μ₀=μ, kwargs...)
    # Set up kernel and streamers with some default values
    streamers = Dict{Playground.Feature, LB.Streaming}(
        Playground.INLET  => NashZeroOrderPressure{T}(T[0, 1], density(p₀)),
        Playground.OUTLET => NashZeroOrderPressure{T}(T[0, 1], convert(T, density(p₀ - Δp))),
        Playground.FLUID  => FluidStreaming(),
        Playground.SOLID  => HalfWayBounceBack()
    )
    result.kernels = Dict(
        Playground.FLUID =>
            FluidKernel(SingleRelaxationTime{T}(1.0/τ), streamers))

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

function lbgk{N, I <: Integer}(lattice::Symbol, dimensions::NTuple{N, I},
                               δt::Number, δx::Number; kwargs...)
    lbgk(getfield(LB, lattice)::LB.Lattice, dimensions, δt, δx; kwargs...)
end
