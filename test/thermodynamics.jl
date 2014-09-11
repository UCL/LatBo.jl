using FactCheck: facts, @fact, roughly
using LatBo: thermodynamics, D2Q9

facts("Verify thermodynamic quantities") do
    context("rho") do
        fᵢ = Float64[1:9...]
        @fact thermodynamics.density(fᵢ) => roughly(10*9/2)
        fᵢ = Float64[1:10...]
        @fact thermodynamics.density(fᵢ) => roughly(10*11/2)
    end

    context("velocity") do
        fᵢ = ones(Float64, 9)
        @fact thermodynamics.velocity(fᵢ, D2Q9.celerities, 1.) .== 0 => all
        fᵢ[1] = 1
        @fact thermodynamics.velocity(fᵢ, D2Q9.celerities, 1.) .== 0 => all
        fᵢ[1] = 2
        actual = thermodynamics.velocity(fᵢ, D2Q9.celerities, 1.)
        @fact actual .== D2Q9.celerities[:, 1] => all
        actual = thermodynamics.velocity(fᵢ, D2Q9.celerities, 2.)
        @fact actual .== (0.5 * D2Q9.celerities[:, 1]) => all
        actual = thermodynamics.velocity(fᵢ, D2Q9.celerities)
        expected = D2Q9.celerities[:, 1] / thermodynamics.density(fᵢ)
        @fact actual .== expected => all
    end

    context("deviatoric tensor") do
        fᵢ = 2*ones(Float64, 9)
        fᵢ⁼ = ones(Float64, 9)
        actual = thermodynamics.deviatoric(fᵢ, fᵢ⁼, D2Q9.celerities, 1.)
        @fact actual .== diagm([3, 3]) => all

        fᵢ[1] = 6
        actual = thermodynamics.deviatoric(fᵢ, fᵢ⁼, D2Q9.celerities, 1.)
        @fact actual .== diagm([3, 3]) => all

        actual = thermodynamics.deviatoric(fᵢ, fᵢ⁼, D2Q9.celerities, 2.)
        @fact actual .== diagm([0, 0]) => all

        actual = thermodynamics.deviatoric(fᵢ, fᵢ⁼, D2Q9.celerities, 0.)
        @fact actual .== diagm([6, 6.]) => all

        fᵢ[2] = 3
        actual = thermodynamics.deviatoric(fᵢ, fᵢ⁼, D2Q9.celerities, 1.)
        c₂ = D2Q9.celerities[:, 2]
        expected = 0.5(c₂ .* transpose(c₂)) + diagm([3, 3])
        @fact actual .== expected => all


    end
end

