module thermodynamics
using LatBo: Lattice

# Density at a given lattice site
density(fᵢ) = sum(fᵢ)
# Momentum at given site
momentum(fᵢ, cᵢ) = vec(sum(cᵢ .* transpose(fᵢ), 2))
# Velocities at a given lattice site
velocity(fᵢ, cᵢ, ρ) = momentum(fᵢ, cᵢ) / ρ
velocity(fᵢ, cᵢ) = velocity(fᵢ, cᵢ, density(fᵢ))

#= Computes the equilibrium particle distributions $f^{eq}$:

    velocity: Macroscopic velocity vector ν at current lattice site
    celerities: d by n matrix of celerities ē for the lattice, with d the dimensionality of the
        lattice and n the number of particle distributions.
    weights: Weights associated with each celerity
    ρ: Density

    $f^{eq}$ = weights .* [ρ + 3ē⋅ν + \frac{9}{2ρ} (ē⋅ν)² - \frac{3}{2ρ} ν⋅ν]$
=#
function equilibrium{T}(ρ::T, velocity::Vector{T}, celerities::Matrix{T}, weights::Vector{T})
    # computes velocity projected on each particle celerity first
    ν_on_ē = celerities.'velocity
    weights .* (
        ρ
        + 3ν_on_ē
        + 9/(2ρ) * (ν_on_ē .* ν_on_ē)
        - 3/(2ρ) * dot(velocity, velocity)
    )
end

function equilibrium{T}(lattice::Lattice, ρ::T, velocity::Vector{T})
  equilibrium(ρ, velocity, lattice.celerities, lattice.weights, lattice.speed_of_sound)
end

end
