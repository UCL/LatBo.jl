using LatBo: Simulation, SandBox
using LatBo.Indices: index
using LatBo.LB: Homogeneous, equilibrium, Initializer
import LatBo.LB.initialize

type MockInitializer <: Initializer
    value
end
function initialize(init::MockInitializer, sim::Simulation, site::Integer)
    sim.populations[:, site] = init.value
end


facts("Initialization") do
    const lattice = :D2Q9
    sim = SandBox(lattice, (4, 4))


    context("Homogeneous momentum and density") do
        const ρ₀ = float64(0.5)
        const μ₀ = Float64[1, 2]
        const feq = equilibrium(:D2Q9, ρ₀, μ₀)
        const start = Int64[3, 3]

        initializer = Homogeneous(ρ₀, μ₀)
        sim.populations[:] = 0
        initialize(initializer, sim, start)
        @fact sim.populations[:, index(sim.indexing, start)] => roughly(feq)
        @fact sum(abs(sim.populations)) => roughly(sum(abs(feq)))
    end

    context("Loop over all sites") do
        sim.playground[:] = 1
        sim.playground[5] = 2
        sim.initializers = {1 => MockInitializer(0), 2 => MockInitializer(1)}

        sim.populations[:] = -1
        initialize(sim)

        @fact any(sim.populations == -1) => false
        @fact sim.populations[:, 5] => roughly(ones(sim.populations[:, 1]))
    end
end
