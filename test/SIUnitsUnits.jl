# Units module using SIUnits.
# SIUnits is currently incompatible with Compose (Keno/SIUnits.jl#19), so this module has been moved
# to tests, where it serves the purpose of an uncommonly long comparison implemention for its more
# primitive replacement.
module Units
using SIUnits: NonSIQuantity, SIQuantity, NonSIUnit, SIUnit, sidims, unit
using SIUnits: quantity, Meter, KiloGram, Pascal, Second

import SIUnits.as
import Base.convert
import Base.isapprox

export LBUnits, as, from, Lenght, Weight, Time, Velocity, Pressure, Viscosity, Density

typealias Length{T}    quantity(T, Meter)
typealias Weight{T}    quantity(T, KiloGram)
typealias Time{T}      quantity(T, Second)
typealias Velocity{T}  quantity(T, Meter/Second)
typealias Pressure{T}  quantity(T, Pascal)
typealias Viscosity{T} quantity(T, Pascal * Second)
typealias Density{T}   quantity(T, KiloGram * Meter^-3)

const mmHg = NonSIUnit{typeof(Pascal), :mmHg}()
convert(::Type{SIQuantity},::typeof(mmHg)) = Pascal / 133.3223684211

# Makes it easier to perform comparison when debugging
function isapprox(x::SIQuantity, y::SIQuantity; kwargs...)
    if unit(x) != unit(y)
        error("Cannot compare objects with different physical dimensions (SIUnits)")
    end
    isapprox(x.val, y.val; kwargs...)
end

# Holds all and everything about units that we know of
immutable LBUnits{T <: Number}
    δt::Time{T}
    δx::Length{T}
    δm::Weight{T}
end
function LBUnits{T <: Number}(::Type{T}, δt::Time, δx::Length, δm::Weight)
    LBUnits{T}(δt, δx, δm)
end
function LBUnits(δt::Time, δx::Length, δm::Weight)
    promoted = promote(δt.val, δx.val, δm.val)
    LBUnits(promoted[1] * unit(δt), promoted[2] * unit(δx), promoted[3] * unit(δm))
end

# Relies on lb density = 1
LBUnits(δt::Time, δx::Length, ρ::Density=1000kg*m^-3) = LBUnits(δt, δx, ρ*δx^3)
function LBUnits{T <: Number}(::Type{T}, δt::Time, δx::Length, ρ::Density=1000kg*m^-3)
    LBUnits(T, δt, δx, ρ*δx^3)
end

# Converts from physical units to (dimensionless) lb units
function as(lb::LBUnits, x::SIUnit)
    if any([u != 0 for u in sidims(x)[4:end]])
        error("Cannot deal with this kind of unit right now")
    end
    dx, dm, dt = sidims(x)
    if dx != 0  lb.δx^(-dx) else 1 end *
        if dm != 0 lb.δm^(-dm) else 1 end *
        if dt != 0 lb.δt^(-dt) else 1 end
end
as(lb::LBUnits, x::SIQuantity) = x * as(lb, unit(x))
function as{T}(lb::LBUnits, x::Array{SIQuantity{T}})
    result = zeros(T, size(x))
    if length(x) > 0
        factor = as(lb, unit(x[1]))
        for i in 1:length(x)
            result[i] = x[i] * factor
        end
    end
    result
end

function from{dx, dm, dt}(lb::LBUnits, ::SIUnit{dx, dm, dt, 0, 0, 0, 0})
    if dx != 0  lb.δx^(dx) else 1 end *
        if dm != 0 lb.δm^(dm) else 1 end *
        if dt != 0 lb.δt^(dt) else 1 end
end
from(lb::LBUnits, x, unit::SIUnit) = x * from(lb, unit)
end
