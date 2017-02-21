# caching of neighbor_index calls with site integer and direction integer
type Cached <: Indexing
    indexing::Indexing
    cache::Matrix
end
function Cached(kern::Indexing, lattice::AbstractLattice)
    const I = typeof(kern).parameters[1]
    result = Cached(kern, zeros(I, (length(lattice), length(kern))))
    for site in 1:length(kern), d in 1:length(lattice)
        result.cache[d, site] = neighbor_index(kern, site, lattice, d)
    end
    result
end
function Cached(kern::Cartesian, lattice::AbstractLattice)
    const I = typeof(kern).parameters[1]
    result = Cached(kern, zeros(I, (length(lattice), length(kern))))
    for site in 1:length(kern), d in 1:length(lattice)
        const coords = gridcoords(kern, site)
        if all(coords .> 1) && all(coords .< size(kern))
            result.cache[d, site] = neighbor_index(kern, site, lattice, d)
        end
    end
    result
end

function cached_cartesian(dims::GridCoords, lattice::AbstractLattice)
    if length(dims) != size(lattice.celerities, 1)
        error("Grid and lattice dimensionality do not match")
    end
    Cached(Cartesian(dims), lattice)
end
function cached_periodic(dims::GridCoords, lattice::AbstractLattice)
    if length(dims) != size(lattice.celerities, 1)
        error("Grid and lattice dimensionality do not match")
    end
    Cached(Periodic(dims), lattice)
end

index(kern::Cached, indices::GridCoords) = index(kern.indexing, indices)
neighbor_index(kern::Cached, site::Integer, ::AbstractLattice, d::Integer) = kern.cache[d, site]
function neighbor_index(kern::Cached, site::GridCoords, ::AbstractLattice, d::Integer)
    kern.cache[d, index(kern.indexing, site)]
end

length(kernel::Cached) = length(kernel.indexing)
size(kernel::Cached) = size(kernel.indexing)
gridcoords(kernel::Cached, index::Integer) = gridcoords(kernel.indexing, index)
function index{I <: Integer}(kernel::Cached, indices::Vararg{I})
    index(kernel.indexing, indices)
end
