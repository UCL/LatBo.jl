module Indices

export index, gridcoords, neighbor_index
import Base: length, size

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

length(d::Indexing) = prod(d.dimensions)
size(d::Indexing) = d.dimensions
function size(d::Indexing, i::Integer)
    @assert i >= 1 and i =< length(d.dimensions)
    d.dimensions[i]
end

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
index(kernel::Indexing, indices::(Integer...)) = index(kernel::Indexing, collect(indices))

function gridcoords(kernel::Indexing, index::Integer)
    result = ones(typeof(index), length(kernel.dimensions))
    index -= 1
    for i in length(kernel.dimensions):-1:1
       result[i] += ifloor(index // kernel.strides[i])
       index %= kernel.strides[i]
    end
    result
end

function neighbor_index(indexing::Indexing, site::GridCoords, direction::GridCoords)
    @assert length(site) == length(direction)
    const neighbor = site + direction
    @assert all(neighbor .> 0) && all(neighbor .<= size(indexing))
    index(indexing, neighbor)
end
function neighbor_index(indexing::Indexing, site::Integer, direction::GridCoords)
    @assert site > 0 && site < length(indexing)
    neighbor_index(indexing, gridcoords(indexing, site), direction)
end

end
