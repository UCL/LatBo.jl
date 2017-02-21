# Adds functionality to deal with transforms from physical to lattice units.
# Mostly, defines two functions:
#
# - `dimensionless` goes from physical to lattice units
# - `physical` goes the other way
#
# Until www.github.com/Keno/SIUnits.jl#19 is resolved, we cannot use SIUnits.jl.
# However, a working implementation does exist in the tests and in the git history.
module Units
export LBunits, dimensionless, physical
# Holds all and everything about units that we know of
immutable LBUnits{T <: Number}
    δt::T
    δx::T
    δm::T
end
function LBUnits{T <: Number}(::Type{T}, δt::Number, δx::Number, δm::Number)
    LBUnits{T}(δt, δx, δm)
end

const DIMENSIONALITY = Dict{Symbol, NTuple{3, Int64}}(
    :time => (1, 0, 0), :length => (0, 1, 0), :weight => (0, 0, 1),
    :velocity => (-1, 1, 0), :pressure => (-2, -1, 1), :viscosity => (-1, -1, 1),
    :density => (0, -3, 1)
)

function dimensionless(lb::LBUnits, dimension::Symbol, x)
    if !haskey(DIMENSIONALITY, dimension)
        error("Don't know how to convert $dimensions to LB units")
    end
    const dt, dx, dm = DIMENSIONALITY[dimension]
    x * (
        if dt != 0 lb.δt^(-dt) else 1 end *
        if dx != 0 lb.δx^(-dx) else 1 end *
        if dm != 0 lb.δm^(-dm) else 1 end
    )
end
function physical(lb::LBUnits, dimension::Symbol, x)
    const dt, dx, dm = DIMENSIONALITY[dimension]
    x * (
        if dt != 0 lb.δt^dt else 1 end *
        if dx != 0 lb.δx^dx else 1 end *
        if dm != 0 lb.δm^dm else 1 end
    )
end
end
