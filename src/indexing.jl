module Indices
# Base type for all indexing kernels
abstract Indexing
typealias GridCoords{I <: Integer} Vector{I}
# No frills indexing
immutable Cartesian{I <: Integer} <: Indexing
    dimensions::GridCoords{I}
    strides::GridCoords{I}
end
Cartesian{I}(dims::GridCoords{I}) = Cartesian{I}(dims, I[strides(zeros(Int8, dims...))...])
# Indices are periodic
immutable Periodic{I <: Integer} <: Indexing
    dimensions::GridCoords{I}
    strides::GridCoords{I}
end
Periodic{I}(dims::GridCoords{I}) = Periodic{I}(dims, I[strides(zeros(Int8, dims...))...])

# A function to retrieve array indices from simulation indices
# Returns a scalar index giving the actual location in space
# Dumps to zero if indices are out of bounds
function index(kernel::Cartesian, indices::GridCoords)
    @assert(size(kernel.dimensions) == size(indices))
    const result = dot(kernel.strides, indices - 1) + 1
    any(indices .< 1) || any(indices .> kernel.dimensions) ? 0: result
end
function index(kernel::Periodic, indices::GridCoords)
    @assert(size(kernel.dimensions) == size(indices))
    dot(kernel.strides, mod(indices - 1, kernel.dimensions)) + 1
end
index(kernel::Indexing, indices::(Integer...)) = index(kernel::Indexing, Integer[indices...])

function gridcoords(kernel::Indexing, index::Integer)
    result = ones(typeof(index), length(kernel.dimensions))
    index -= 1
    for i in length(kernel.dimensions):-1:1
       result[i] += ifloor(index // kernel.strides[i])
       index %= kernel.strides[i]
    end
    result
end

end
