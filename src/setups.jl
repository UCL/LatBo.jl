export lbgkL
using SIUnits

typealias Length{T} SIUnits.SIQuantity{T,1,0,0,0,0,0,0}
typealias Weight{T} SIUnits.SIQuantity{T,0,1,0,0,0,0,0}
typealias Time{T} SIUnits.SIQuantity{T,0,0,1,0,0,0,0}
typealias Velocity{T} SIUnits.SIQuantity{T,1,0,-1,0,0,0,0}
typealias Pressure{T} SIUnits.SIQuantity{T,-1,1,-2,0,0,0,0}
typealias Viscosity{T} SIUnits.SIQuantity{T,-1,1,-3,0,0,0,0}
typealias Density{T} SIUnits.SIQuantity{T,-3,1,0,0,0,0,0}

const mmHg = SIUnits.NonSIUnit{typeof(SIUnits.Pascal), :mmHg}()
convert(::Type{SIUnits.SIQuantity}, ::typeof(mmHg)) = SIUnits.Pascal / 133.3223684211

# Sets up calculation with standard lbgk
function lbgk{T <: FloatingPoint, I <| FloatingPoint}(lattice::LB.Lattice, dimensions::(Integer...),
        δt::Time, viscosity::Viscosity;
        ρ₀::Density=1e3kg/m^3, μ₀::Velocity=[0, 0]m/s, p₀::Pressure=80.0mmHg, kwargs...)

    const μₗ = μ₀δt/δx
    const pₗ = p₀/ρ₀*δx^2/δt
    const cₛ² = 1./3.
    const cₛ = sqrt(cₛ²)

    density(p) = 1. + convert()

    result = SandBox(lattice, dimensions; ρ₀=ρ₀/ρₗ, μ₀=μ₀/μₗ, kwargs...)
    result.kernels = {Playground.FLUID => SingleRelaxationTime{T, V}(1.0/τ)}
end

