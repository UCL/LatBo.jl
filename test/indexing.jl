using LatBo.Indices: Cartesian, Periodic, index

facts("Indexing") do

    context("Cartesian") do
        @fact index(Cartesian([10, 10, 10]), [1, 2, 3]) => [1, 2, 3]
        @fact index(Cartesian([10, 10, 30]), [8, 10, 20]) => [8, 10, 20]
        @fact index(Cartesian([10, 10, 30]), [0, 10, 20]) => [0, 0, 0]
        @fact index(Cartesian([10, 10, 30]), [10, 11, 20]) => [0, 0, 0]
    end

    context("Periodic") do
        @fact index(Periodic([64, 64, 64]), [8, 16, 32]) => [8, 16, 32]
        @fact index(Periodic([5, 64, 64]), [8, 16, 32]) => [3, 16, 32]
        # Between 1 and n, not 0 and n-1. Julia indices start at 1
        @fact index(Periodic([64, 16, 64]), [8, 0, 32]) => [8, 16, 32]
        @fact index(Periodic([64, 64, 32]), [8, 16, 32]) => [8, 16, 32]
        @fact index(Periodic([64, 64, 32]), [8, 16, 33]) => [8, 16, 1]
    end
end
