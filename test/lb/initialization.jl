using LatBo.LatticeBoltzmann: Homogeneous, initialize, equilibrium

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
        @fact sim.populations[:, start...] => roughly(feq)
        @fact sum(abs(sim.populations)) => roughly(sum(abs(feq)))
    end
end
