using LatBo: Simulation
using LatBo.LB: D2Q9, D3Q19, Lattice
using LatBo.Indices: Indexing

facts("Lattice direction inversion") do
    for lattice_name in [:D2Q9, :D3Q19]
        context("for $(lattice_name)") do
            lattice = eval(lattice_name)
            ncomponents = size(lattice.celerities, 2)

            @fact lattice.inversion .>= 1 --> all
            @fact lattice.inversion .<= ncomponents --> all
            @fact length(unique(lattice.inversion)) --> ncomponents
            @fact lattice.celerities[:, lattice.inversion] --> -lattice.celerities
        end
    end
end

type FakeSimulation <: Simulation
    indexing::Indexing
    lattice::Lattice
end

facts("Compute neighbor index from sim/lattice object") do
    sim = FakeSimulation(Cartesian([10, 7]), D2Q9)
    for from = ([2, 3], [5, 6]), direction = 1:size(sim.lattice.celerities, 2)
        site = index(sim.indexing, from)
        to = neighbor_index(sim, site, direction)
        @fact gridcoords(sim.indexing, to) - from --> sim.lattice.celerities[:, direction]
    end
    @fact_throws AssertionError neighbor_index(sim, 0, 1)
    @fact_throws AssertionError neighbor_index(sim, length(sim.indexing) + 1, 1)
end
