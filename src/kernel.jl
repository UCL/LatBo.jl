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
    dimensions::Vector{Int64}
end
# Indices are periodic
immutable Periodic <: Indexing
    dimensions::Vector{Int64}
end

# A function to retrieve array indices from simulation indices
# Dumps to zero if indices are out of bounds
function index{T1 <: Int, T2 <: Int}(dimensions::Vector{T1}, indices::Vector{T2})
    @assert(size(dimensions) == size(indices))
    any(indices .< 1) || any(indices .> dimensions) ?
        zeros(T2, size(indices)): indices
end
index(kernel::Indexing, indices) = index(kernel.dimensions, indices)
index{T}(kernel::Periodic, indices::Array{T}) = 1 + mod(indices - 1, kernel.dimensions)

end
