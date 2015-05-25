# Indices are periodic
immutable Periodic{I <: Integer} <: Indexing
    dimensions::GridCoords{I}
    strides::GridCoords{I}
end
Periodic{I}(dims::GridCoords{I}) = Periodic{I}(dims, I[strides(zeros(Int8, dims...))...])

function index(kernel::Periodic, indices::GridCoords)
    @assert(size(kernel.dimensions) == size(indices))
    dot(kernel.strides, mod(indices - 1, kernel.dimensions)) + 1
end
function neighbor_index(indexing::Periodic, site::GridCoords, direction::GridCoords)
    @assert length(site) == length(direction)
    index(indexing, site + direction)
end
