using LatBo.Units: LBUnits, as, from
using SIUnits: SIQuantity, SIUnit, unit
using SIUnits.ShortUnits

facts("Physical to and from LB units") do
    dimensionless{T<:Number}(::Type{T}) = !(T <: SIQuantity) && !(T <: SIUnit)
    dimensionless(x::Number) = dimensionless(typeof(x))
    dimensionless(x::Array) = dimensionless(eltype(x))

    δt, δx, δm = 0.1, 2., 30.
    lbunits = LBUnits{Float64}(δt * s, δx * m, δm * kg)

    context("P to LB") do
        @fact as(lbunits, 0.1m) => dimensionless
        @fact as(lbunits, 2m) => roughly(2 / δx)
        @fact as(lbunits, 1m) => roughly(1 / δx)
        @fact as(lbunits, 1/m) => roughly(δx)
        @fact as(lbunits, 1s) => dimensionless
        @fact as(lbunits, 1s) => roughly(1 / δt)
        @fact as(lbunits, 1m/s) => roughly(δt/δx)
        @fact as(lbunits, 0.5kg/m^3/s) => roughly(0.5δt*δx^3/δm)
        @fact as(lbunits, [0.5, 0.6]*(kg*m^-3/s)) => roughly([0.5, 0.6]δt*δx^3/δm)
        @fact as(lbunits, [0.5, 0.6]*(kg*m^-3/s)) => dimensionless
        @fact_throws as(lbunits, 0.1Pascal * mol)
    end

    context("LB to P") do
        @fact from(lbunits, m) => roughly(lbunits.δx)
        @fact unit(from(lbunits, 2, m)) => m
        @fact from(lbunits, 2, m) => roughly(2.lbunits.δx)
        some_dim = lbunits.δm * lbunits.δx^-3 * lbunits.δt^-1
        @fact from(lbunits, 0.5, kg/m^3/s) => roughly(0.5 * some_dim)
        @fact from(lbunits, [0.5, 0.6], kg*m^-3/s) => roughly([0.5, 0.6] * some_dim)
    end

    context("P to LB to P") do
        for value in [1m, 2m, 3kg/s/m^2]
            @fact from(lbunits, as(lbunits, value), unit(value)) => roughly(value)
        end
    end
end
