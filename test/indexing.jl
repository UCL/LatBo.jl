using LatBo.Indices: Cartesian, Periodic, index, gridcoords, neighbor_index

facts("Grid coordinates to array index") do
    context("Cartesian") do
        @fact index(Cartesian([10, 10, 10]), [1, 2, 3]) => 0 + 10 * 1 + 100 * 2 + 1
        @fact index(Cartesian([10, 10, 30]), [8, 10, 20]) => 7 + 10 * 9 + 100 * 19 + 1
        @fact index(Cartesian([10, 10, 30]), [0, 10, 20]) => 0
        @fact index(Cartesian([10, 10, 30]), [10, 11, 20]) => 0
    end

    context("Periodic") do
        @fact index(Periodic([64, 64, 64]), [8, 16, 32]) => 7 + 64 * 15 + 64 * 64 * 31 + 1
        @fact index(Periodic([5, 64, 64]), [8, 16, 32]) => 2 + 5 * 15 + 5 * 64 * 31 + 1
        # Between 1 and n, not 0 and n-1. Julia indices start at 1
        @fact index(Periodic([64, 16, 64]), [8, 0, 32]) => dot([7, 15, 31], [1, 64, 16 * 64]) + 1
        @fact index(Periodic([64, 64, 32]), [8, 16, 32]) => dot([7, 15, 31], [1, 64, 64 * 64]) + 1
        @fact index(Periodic([64, 64, 32]), [8, 16, 33]) => dot([7, 15, 0], [1, 64, 64 * 64]) + 1
    end
end

facts("Array index to grid coordinates") do
    @fact gridcoords(Cartesian([10, 20, 30]), 1) => [1, 1, 1]
    @fact gridcoords(Cartesian([10, 20, 30]), 10 * 20 * 30) => [10, 20, 30]
    @fact gridcoords(Cartesian([10, 20, 30]), 1 + 4 + 10 * 4 + 200 * 4) => [5, 5, 5]
end

facts("Compute neighbor index") do
    cart = Cartesian([10, 7])
    for from = ([2, 3], [5, 6]), direction = ([1, 0], [0, 1], [1, -1], [0, 0])
        site = index(cart, from)
        to = neighbor_index(cart, site, direction)
        @fact gridcoords(cart, to) - from => direction
    end
    @fact_throws ErrorException neighbor_index(cart, 0, [0, 0])
    @fact_throws ErrorException neighbor_index(cart, length(cart) + 1, [0, 0])
end
