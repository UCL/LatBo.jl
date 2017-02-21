type LocalQuantities{T <: AbstractFloat} <: AbstractLocalQuantities
    density::T
    momentum::DenseVector{T}
    velocity::DenseVector{T}
    feq::DenseVector{T}
end

function LocalQuantities{T}(fᵢ::DenseVector{T}, lattice::Lattice{T})
    const ρ = density(fᵢ)
    const μ = momentum(fᵢ, lattice.celerities)
    const ν = velocity(μ, ρ)
    const feq = equilibrium(ρ, μ, lattice.celerities, lattice.weights)
    LocalQuantities{T}(ρ, μ, ν, feq)
end

# Pre-allocates a local quantity aggregator
function LocalQuantities(lattice::Lattice)
    const T = eltype(lattice.weights)

    LocalQuantities{T}(
        convert(T, 0), zeros(T, ndims(lattice)),
        zeros(T, ndims(lattice)), zeros(T, length(lattice))
    )
end

# Recomputes local-quantities in-place
function LocalQuantities!{T}(out::LocalQuantities{T}, fᵢ::DenseVector{T}, lattice::Lattice{T})
    out.density = density(fᵢ)
    momentum!(out.momentum, fᵢ, lattice.celerities)
    velocity!(out.velocity, out.momentum, out.density)
    equilibrium!(out.feq, out.density, out.momentum, lattice.celerities, lattice.weights)
    out
end
