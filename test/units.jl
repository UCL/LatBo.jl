using SIUnits: SIQuantity, SIUnit, unit
using SIUnits.ShortUnits
using LatBo

include("SIUnitsUnits.jl")
using LatBoTests.Units: as, from

facts("Physical to and from LB units using SIUnits") do
    is_dimensionless{T<:Number}(::Type{T}) = !(T <: SIQuantity) && !(T <: SIUnit)
    is_dimensionless(x::Number) = is_dimensionless(typeof(x))
    is_dimensionless(x::Array) = is_dimensionless(eltype(x))

    δt, δx, δm = 0.1, 2., 30.
    siunits = LatBoTests.Units.LBUnits{Float64}(δt * s, δx * m, δm * kg)

    context("P to LB") do
        @fact as(siunits, 0.1m) --> is_dimensionless
        @fact as(siunits, 2m) --> roughly(2 / δx)
        @fact as(siunits, 1m) --> roughly(1 / δx)
        @fact as(siunits, 1/m) --> roughly(δx)
        @fact as(siunits, 1s) --> is_dimensionless
        @fact as(siunits, 1s) --> roughly(1 / δt)
        @fact as(siunits, 1m/s) --> roughly(δt/δx)
        @fact as(siunits, 0.5kg/m^3/s) --> roughly(0.5δt*δx^3/δm)
        @fact as(siunits, [0.5, 0.6]*(kg*m^-3/s)) --> roughly([0.5, 0.6]δt*δx^3/δm)
        @fact as(siunits, [0.5, 0.6]*(kg*m^-3/s)) --> is_dimensionless
        @fact_throws as(siunits, 0.1Pascal * mol)
    end

    context("LB to P") do
        @fact from(siunits, m) --> roughly(siunits.δx)
        @fact unit(from(siunits, 2, m)) --> m
        @fact from(siunits, 2, m) --> roughly(2.siunits.δx)
        some_dim = siunits.δm * siunits.δx^-3 * siunits.δt^-1
        @fact from(siunits, 0.5, kg/m^3/s) --> roughly(0.5 * some_dim)
        @fact from(siunits, [0.5, 0.6], kg*m^-3/s) --> roughly([0.5, 0.6] * some_dim)
    end

    context("P to LB to P") do
        for value in [1m, 2m, 3kg/s/m^2]
            @fact from(siunits, as(siunits, value), unit(value)) --> roughly(value)
        end
    end
end


facts("Physical to and from LB units") do
    δt, δx, δm = 0.1, 2., 30.
    lbunits = LatBo.Units.LBUnits{Float64}(δt, δx, δm)
    siunits = LatBoTests.Units.LBUnits{Float64}(δt * s, δx * m, δm * kg)

    equivs = [
        (:time, s), (:length, m), (:weight, kg), (:density, 1kg*m^-3),
        (:velocity, 1m/s), (:pressure, Pa), (:viscosity, 1Pa * s)
    ]
    goaround(dim::Symbol, x) = LatBo.Units.physical(
        lbunits, dim, LatBo.Units.dimensionless(lbunits, dim, x))
    for (lb, si) in equivs
        @fact LatBo.Units.dimensionless(lbunits, lb, 2) --> roughly(as(siunits, 2 * si))
        r = 0.5 * randn(10, 5) + 0.5
        @fact LatBo.Units.dimensionless(lbunits, lb, r) --> roughly(as(siunits, r * si))

        @fact goaround(lb, 2) --> roughly(2)
        @fact goaround(lb, r) --> roughly(r)
    end

    @fact_throws LatBo.Units.dimensionles(lbunits, :notaunit, 2)
end
