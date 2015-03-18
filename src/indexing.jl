module Indices
# Base type for all indexing kernels
abstract Indexing
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
index(kernel::Indexing, indices) = indices
index(kernel::Cartesian, indices) = index(kernel.dimensions, indices)
index{T <: Int}(kernel::Periodic, indices::Array{T}) = 1 + mod(indices - 1, kernel.dimensions)
end
