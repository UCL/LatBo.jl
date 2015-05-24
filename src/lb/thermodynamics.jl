# Density at a given lattice site
density(fᵢ::DenseVector) = sum(fᵢ)
function density(fᵢ::DenseArray)
    shape = size(fᵢ)[2:end]
    reshape([sum(fᵢ[:, i]) for i = 1:prod(shape)], shape)
end
density(sim::Simulation) = density(sim.populations)
# Momentum at given site
momentum(fᵢ::DenseVector, cᵢ::DenseMatrix) = cᵢ * fᵢ
# Velocities at a given lattice site
velocity(μ::DenseVector, ρ::Number) = μ / ρ
velocity(fᵢ::DenseVector, cᵢ::DenseMatrix, ρ::Number) = velocity(momentum(fᵢ, cᵢ), ρ)
velocity(fᵢ::DenseVector, cᵢ::DenseMatrix) = velocity(fᵢ, cᵢ, density(fᵢ))

# Momentum and velocity functions that take the full grid, or a simulation object
for name in [:momentum, :velocity]
    @eval begin
        function $name(fᵢ::DenseArray, cᵢ::DenseMatrix)
            const shape = tuple(size(cᵢ, 1), size(fᵢ)[2:end]...)
            result = zeros(eltype(fᵢ), shape)
            for i = 1:prod(shape[2:end])
                result[1:shape[1], i] = $name(fᵢ[:, i], cᵢ)
            end
            result
        end
        $name(sim::Simulation) = $name(sim.populations, sim.lattice.celerities)
    end
end

#= Computes the equilibrium particle distributions $f^{eq}$:

    momentum: Macroscopic momentum μ at current lattice site
    celerities: d by n matrix of celerities ē for the lattice, with d the dimensionality of the
        lattice and n the number of particle distributions.
    weights: Weights associated with each celerity
    ρ: Density

    $f^{eq}$ = weights .* [ρ + 3ē⋅μ + \frac{9}{2ρ} (ē⋅μ)² - \frac{3}{2ρ} μ⋅μ]$
=#
function equilibrium{T, I}(
        ρ::T, momentum::DenseVector{T}, celerities::DenseMatrix{I}, weights::DenseVector{T})
    # computes momentum projected on each particle celerity first
    @assert length(momentum) == size(celerities, 1)
    @assert length(weights) == size(celerities, 2)
    μ_on_ē = celerities.'momentum
    weights .* (
        ρ
        + 3μ_on_ē
        + 9/(2ρ) * (μ_on_ē .* μ_on_ē)
        - 3/(2ρ) * dot(momentum, momentum)
    )
end
equilibrium{T, I}(lattice::Lattice{T, I}, fᵢ::DenseVector{T}) =
    equilibrium(density(fᵢ), momentum(fᵢ, lattice.celerities), lattice.celerities, lattice.weights)

equilibrium{T}(lattice::Lattice, ρ::T, momentum::DenseVector{T}) =
    equilibrium(ρ, momentum, lattice.celerities, lattice.weights)
equilibrium{T}(lattice::Symbol, ρ::T, momentum::DenseVector{T}) =
    equilibrium(getfield(LB, lattice), ρ, momentum)

immutable type LocalQuantities{T <: FloatingPoint, I <: Int}
    density::T
    momentum::DenseVector{T}
    velocity::DenseVector{T}
    feq::DenseVector{T}

    function LocalQuantities(fᵢ::DenseVector{T}, lattice::Lattice{T, I})
        const ρ = density(fᵢ)
        const μ = momentum(fᵢ, lattice.celerities)
        const ν = velocity(μ, ρ)
        const feq = equilibrium(ρ, μ, lattice.celerities, lattice.weights)
        new(ρ, μ, ν, feq)
    end
end

# Convenience calls to specify types as arguments
LocalQuantities{T <: FloatingPoint, I <: Integer}(::Type{T}, ::Type{I}, args...) =
    LocalQuantities{T, I}(args...)
function LocalQuantities(types::(Type, Type), args...)
    @assert types[1] <: FloatingPoint
    @assert types[2] <: Integer
    LocalQuantities(types..., args...)
end
