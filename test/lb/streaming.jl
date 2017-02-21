using LatBo.LB: HalfWayBounceBack, FluidStreaming, streaming
using LatBo: SandBox
using LatBo.Indices: index
facts("Streaming kernels") do
    sim = SandBox(:D2Q9, (4, 4))

    context("Fluid to Fluid streaming") do
        sim.populations[:] = 0
        sim.next_populations[:] = 0
        const cᵢ = sim.lattice.celerities
        const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
        const start = index(sim.indexing, (3, 3))
        const finish = index(sim.indexing, (2, 4))
        sim.populations[direction, start] = 1
        @fact sim.next_populations .== 0 --> all
        # Perform streaming
        streaming(FluidStreaming(), sim, start, finish, direction)
        # Check it went to the right place and only there
        @fact sim.next_populations[direction, finish] --> roughly(1)
        @fact sum(abs(sim.next_populations)) --> roughly(1)
        # Check populations was not affected
        @fact sim.populations[direction, start] --> roughly(1)
        @fact sum(abs(sim.populations)) --> roughly(1)
    end

    context("Halfway bounce-back streaming") do
        sim.populations[:] = 0
        sim.next_populations[:] = 0
        const cᵢ = sim.lattice.celerities
        const direction = find([all(cᵢ[:, i] .== [-1, 1]) for i in 1:size(cᵢ, 2)])[1]
        const invdir = find([all(cᵢ[:, i] .== [1, -1]) for i in 1:size(cᵢ, 2)])[1]
        const start = index(sim.indexing, (3, 3))

        @fact sim.lattice.inversion[direction] --> invdir
        sim.populations[direction, start...] = 1
        @fact sim.next_populations .== 0 --> all
        # Perform streaming
        streaming(HalfWayBounceBack(), sim, start, direction)
        # Check it went to the right place and only there
        @fact sim.next_populations[invdir, start...] --> roughly(1)
        @fact sum(abs(sim.next_populations)) --> roughly(1)
        @fact sim.populations[direction, start...] --> roughly(1)
        @fact sum(abs(sim.populations)) --> roughly(1)
    end
end
