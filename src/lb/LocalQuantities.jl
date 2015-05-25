type LocalQuantities{T <: FloatingPoint, I <: Int} <: AbstractLocalQuantities
    density::T
    momentum::DenseVector{T}
    velocity::DenseVector{T}
    feq::DenseVector{T}
end

function LocalQuantities{T, I}(fᵢ::DenseVector{T}, lattice::Lattice{T, I})
    const ρ = density(fᵢ)
    const μ = momentum(fᵢ, lattice.celerities)
    const ν = velocity(μ, ρ)
    const feq = equilibrium(ρ, μ, lattice.celerities, lattice.weights)
    LocalQuantities{T, I}(ρ, μ, ν, feq)
end

# Pre-allocates a local quantity aggregator
function LocalQuantities(lattice::Lattice)
    const T = eltype(lattice.weights)
    const I = eltype(lattice.inversion)

    LocalQuantities{T, I}(
        convert(T, 0), zeros(T, ndims(lattice)),
        zeros(T, ndims(lattice)), zeros(T, length(lattice))
    )
end

# Convenience calls to specify types as arguments
LocalQuantities{T <: FloatingPoint, I <: Integer}(::Type{T}, ::Type{I}, args...) =
    LocalQuantities(args...)::LocalQuantities{T, I}
function LocalQuantities(types::(Type, Type), args...)
    @assert types[1] <: FloatingPoint
    @assert types[2] <: Integer
    LocalQuantities(types..., args...)
end

# Recomputes local-quantities in-place
function LocalQuantities!{T, I}(
        out::LocalQuantities{T, I}, fᵢ::DenseVector{T}, lattice::Lattice{T, I})
    out.density = density(fᵢ)
    momentum!(out.momentum, fᵢ, lattice.celerities)
    velocity!(out.velocity, out.momentum, out.density)
    equilibrium!(out.feq, out.density, out.momentum, lattice.celerities, lattice.weights)
    out
end

