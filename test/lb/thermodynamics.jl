using LatBo.LB: velocity, momentum, density, equilibrium, D2Q9, D3Q19, LocalQuantities

facts("Thermodynamic quantities and functions") do
    for lattice_name in [:D2Q9, :D3Q19]
        context("for $(lattice_name)") do
            lattice = eval(lattice_name)
            ncomponents = size(lattice.celerities, 2)

            context("rho is a geometric series") do
                @fact density(Float64[1:9...]) => roughly(10*9/2)
                @fact density(Float64[1:10...]) => roughly(10*11/2)
            end

            context("homogeneous populations sum to zero momentum") do
                fᵢ = ones(Float64, ncomponents)
                @fact momentum(fᵢ, lattice.celerities) .== 0 => all
                fᵢ[1] = 3 # This is the zero momentum component
                @fact momentum(fᵢ, lattice.celerities) .== 0 => all
                fᵢ[2] = 2 # Now for a negative test
                @fact momentum(fᵢ, lattice.celerities) .!= 0 => any
            end

            context("velocity from momentum and density ") do
                fᵢ = convert(Array{Float64}, 1.0 + rand(ncomponents))
                actual = velocity(fᵢ, lattice.celerities, 1.)
                half = velocity(fᵢ, lattice.celerities, 0.5)
                μ = momentum(fᵢ, lattice.celerities)
                @fact actual => roughly(half * 0.5)
                @fact actual => roughly(μ)
            end

            context("equilibrium function") do
                wᵢ = Float64[1, 2, 3]
                ē = Int64[1 0 1; 0 1 1]
                μ = Float64[1, 2]
                ρ = 1.1

                @fact equilibrium(ρ, zeros(2), ē, wᵢ) => roughly(wᵢ * ρ)
                @fact equilibrium(ρ, μ, zeros(Int64, 2, 3), wᵢ) => roughly(wᵢ * (ρ - 1.5  * 5./ρ))

                fᵉ = wᵢ .* (ρ - 1.5dot(μ, μ)/ρ + 4.5(ē.'μ).^2/ρ + 3ē.'μ)
                @fact equilibrium(ρ, μ, ē, wᵢ) => roughly(fᵉ)
            end

            context("local quantities aggregator") do
               from = Int64[1, 2, 3]
               fᵢ = convert(Array{Float64}, 1.0 + rand(Float64, ncomponents))
               μ = momentum(fᵢ, lattice.celerities)
               ν = velocity(fᵢ, lattice.celerities)
               ρ = density(fᵢ)
               feq = equilibrium(lattice, ρ, μ)

               quants = LocalQuantities{Float64, Int64}(from, fᵢ, lattice)
               @fact quants.from => from
               @fact quants.density => roughly(ρ)
               @fact quants.momentum => roughly(μ)
               @fact quants.velocity => roughly(ν)
               @fact quants.feq => roughly(feq)
            end
        end
    end
end
