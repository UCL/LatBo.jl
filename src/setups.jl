export lbgkL
using SIUnits

typealias Length{T} SIUnits.SIQuantity{T,1,0,0,0,0,0,0}
typealias Weight{T} SIUnits.SIQuantity{T,0,1,0,0,0,0,0}
typealias Time{T} SIUnits.SIQuantity{T,0,0,1,0,0,0,0}
typealias Velocity{T} SIUnits.SIQuantity{T,1,0,-1,0,0,0,0}
typealias Pressure{T} SIUnits.SIQuantity{T,-1,1,-2,0,0,0,0}
typealias Viscosity{T} SIUnits.SIQuantity{T,-1,1,-3,0,0,0,0}
typealias Density{T} SIUnits.SIQuantity{T,-3,1,0,0,0,0,0}

# Sets up calculation with standard lbgk
function lbgk(lattice::LB.Lattice, dimensions::(Integer...),
        δt::Time, viscosity::Viscosity; ρ₀::Density=1e3kg/m^3, μ₀::Velocity=[0, 0]m/s, kwargs...)
    result = SandBox(lattice, dimensions; kwargs)
    result.kernels
end

