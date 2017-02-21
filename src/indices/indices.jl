module Indices
using ..AbstractLattice

export index, gridcoords, neighbor_index
import Base: length, size

# Base type for all indexing kernels
abstract Indexing
typealias GridCoords{I <: Integer} Vector{I}

include("cartesian.jl")
include("periodic.jl")

length(d::Indexing) = prod(d.dimensions)
size(d::Indexing) = d.dimensions
function size(d::Indexing, i::Integer)
    @assert i ≥ 1 && i ≤ length(d.dimensions)
    d.dimensions[i]
end

function index{N, I <: Integer}(kernel::Indexing, indices::NTuple{N, I})
    index(kernel::Indexing, collect(indices))
end

function gridcoords(kernel::Indexing, index::Integer)
    result = ones(typeof(index), length(kernel.dimensions))
    index -= 1
    for i in length(kernel.dimensions):-1:1
       result[i] += floor(typeof(index), index // kernel.strides[i])
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

include("cached.jl")

end
