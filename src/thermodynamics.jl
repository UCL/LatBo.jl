module thermodynamics
using ..Lattice

# Density at a given lattice site
density(fᵢ::Vector) = sum(fᵢ)
# Momentum at given site
momentum(fᵢ::Vector, cᵢ::Matrix) = vec(sum(cᵢ .* transpose(fᵢ), 2))
# Velocities at a given lattice site
velocity(μ::Vector, ρ::Number) = μ / ρ
velocity(fᵢ::Vector, cᵢ::Matrix, ρ::Number) = velocity(momentum(fᵢ, cᵢ), ρ)
velocity(fᵢ::Vector, cᵢ::Matrix) = velocity(fᵢ, cᵢ, density(fᵢ))

#= Computes the equilibrium particle distributions $f^{eq}$:

    momentum: Macroscopic momentum μ at current lattice site
    celerities: d by n matrix of celerities ē for the lattice, with d the dimensionality of the
        lattice and n the number of particle distributions.
    weights: Weights associated with each celerity
    ρ: Density

    $f^{eq}$ = weights .* [ρ + 3ē⋅μ + \frac{9}{2ρ} (ē⋅μ)² - \frac{3}{2ρ} μ⋅μ]$
=#
function equilibrium{T, I}(ρ::T, momentum::Vector{T}, celerities::Matrix{I}, weights::Vector{T})
    # computes momentum projected on each particle celerity first
    μ_on_ē = celerities.'momentum
    weights .* (
        ρ
        + 3μ_on_ē
        + 9/(2ρ) * (μ_on_ē .* μ_on_ē)
        - 3/(2ρ) * dot(momentum, momentum)
    )
end

equilibrium{T}(lattice::Lattice, ρ::T, momentum::Vector{T}) =
    equilibrium(ρ, momentum, lattice.celerities, lattice.weights)

immutable type LocalQuantities{T <: Real, I <: Int}
    from::Vector{I}
    density::T
    momentum::Vector{T}
    velocity::Vector{T}
    feq::Vector{T}

    function LocalQuantities(from::Vector{I}, fᵢ::Vector{T}, lattice::Lattice{T, I})
        const ρ = thermodynamics.density(fᵢ)
        const μ = thermodynamics.momentum(fᵢ, lattice.celerities)
        const ν = thermodynamics.velocity(μ, ρ)
        const feq = thermodynamics.equilibrium(ρ, μ, lattice.celerities, lattice.weights)
        new(from, ρ, μ, ν, feq)
    end
end


end
