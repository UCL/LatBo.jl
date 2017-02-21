using LatBo.Indices.neighbor_index
facts("Checking caching") do
    lattice = LatBo.LB.D2Q9
    context("Input throws if incorrect") do
        @fact_throws ErrorException LatBo.Indices.cached_cartesian([5, 4, 2], lattice)
        @fact_throws ErrorException LatBo.Indices.cached_periodic([5, 4, 2], lattice)
    end

    context("Cached Cartesian") do
        cartesian = LatBo.Indices.Cartesian([5, 4])
        cached = LatBo.Indices.Cached(cartesian, lattice)

        @fact size(cached.cache) --> (length(lattice), length(cartesian))
        @fact typeof(cached.indexing) --> typeof(cartesian)
        @fact size(cached) --> size(cartesian)
        @fact length(cached) --> length(cartesian)
        for i = 1:length(cartesian), d = 1:length(lattice)
            const coords = gridcoords(cartesian, i)
            if all(coords .> 1) && all(coords .< size(cartesian))
                const expected = neighbor_index(cartesian, i, lattice, d)
                @fact neighbor_index(cached, i, lattice, d) --> expected
            else
                @fact neighbor_index(cached, i, lattice, d) --> 0
            end
        end
    end

    context("Cached Periodic") do
        periodic = LatBo.Indices.Periodic([5, 4])
        cached = LatBo.Indices.Cached(periodic, lattice)
        @fact size(cached.cache) --> (length(lattice), length(periodic))
        @fact typeof(cached.indexing) --> typeof(periodic)
        @fact size(cached) --> size(periodic)
        @fact length(cached) --> length(periodic)
        for i = 1:length(periodic), d = 1:length(lattice)
            @fact neighbor_index(cached, i, lattice, d) --> neighbor_index(periodic, i, lattice, d)
        end
    end
end
