# No frills indexing
immutable Cartesian{I <: Integer} <: Indexing
    dimensions::GridCoords{I}
    strides::GridCoords{I}
end
Cartesian{I}(dims::GridCoords{I}) = Cartesian{I}(dims, I[strides(zeros(Int8, dims...))...])

# A function to retrieve array indices from simulation indices
# Returns a scalar index giving the actual location in space
# Dumps to zero if indices are out of bounds
function index(kernel::Cartesian, indices::GridCoords)
    @assert(size(kernel.dimensions) == size(indices))
    const result = dot(kernel.strides, indices - 1) + 1
    any(indices .< 1) || any(indices .> kernel.dimensions) ? 0: result
end
