using Base: Cartesian
module kernel

# Base type for all kernally stuff
abstract Kernel
# Base type for all collision kernels
abstract Collision <: Kernel
type SingleRelaxationTime{T <: Real} <: Collision
    # Inverse of the relaxation time
    τ⁻¹::T
end

# Defines collision kernels for SRT
collision{T}(τ⁻¹::T, fᵢ::Vector{T}, feq::Vector{T}) = τ⁻¹ * (feq - fᵢ)
collision{T}(k::SingleRelaxationTime{T}, fᵢ::Vector{T}, feq::Vector{T}) = collision(k.τ⁻¹, feq, fᵢ)

# Base type for all indexing kernels
abstract Indexing <: Kernel
# No frills indexing
immutable Cartesian <: Indexing
end
# Indices are periodic
immutable Periodic <: Indexing
    dimensions::(Int64...)
end

# A function to retrieve array indices from simulation indices
index(kernel::Indexing, indices) = indices
function index{T}(kernel::Periodic, indices::Array{T})
  T[1 + mod(i - 1, convert(T, d)) for (i, d) in zip(indices, kernel.dimensions)]
end

end
