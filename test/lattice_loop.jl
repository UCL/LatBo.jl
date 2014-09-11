using FactCheck: facts, @fact, roughly, context
using LatBo: SingleRelaxationTime, lattice_loop

facts("Single relaxation time initialization") do
    context("2D") do
        sim = SingleRelaxationTime(0.5, (5, 5))
        sim.populations = rand(size(sim.populations)...)
        populations = copy(sim.populations)
        sim.playground = rand(size(sim.playground)...)
        lattice_loop(sim) do indices, fᵢ, feature
            @fact fᵢ => roughly(sim.populations[:, indices...])
            @fact feature => roughly(sim.playground[indices...])
            populations[:, indices...] = 0
        end
        @fact populations => roughly(zeros(size(populations)...))
    end
    context("3D") do
        sim = SingleRelaxationTime(0.5, (3, 3, 3))
        sim.populations = rand(size(sim.populations)...)
        populations = copy(sim.populations)
        sim.playground = rand(size(sim.playground)...)
        lattice_loop(sim) do indices, fᵢ, feature
            @fact fᵢ => roughly(sim.populations[:, indices...])
            @fact feature => roughly(sim.playground[indices...])
            populations[:, indices...] = 0
        end
        @fact populations => roughly(zeros(size(populations)...))
    end
end
