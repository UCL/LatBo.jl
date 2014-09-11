using FactCheck: facts, @fact, roughly, @fact_throws, context, not
using LatBo: SingleRelaxationTime, run_lb
using LatBo.playground: FLUID, INLET, SOLID

facts("Goes one step through main loop") do
    sim = SingleRelaxationTime(1., (6, 6))

    context("no boundaries, no inlets") do
        sim.playground[:] = FLUID
        sim.populations[:] = 0
        sim.populations[1, :, :] = 1
        sim.populations[2, 1:2:end, 1:2:end] = 1
        sim.populations[1, 1:2:end, 1:2:end] = 0

        # Run one step
        run_lb(sim)

        # Check things streaming right
        pop = sim.populations[2, 2:2:end, 2:2:end]
        @fact pop - pop[1] => roughly(zeros(size(pop)...))
        pop = sim.populations[2, 1:2:end, 1:end]
        @fact pop - pop[1] => roughly(zeros(size(pop)...))
        @fact sim.populations[2, 2, 1] => not(roughly(sim.populations[2, 2, 2]))

        # Check that zero velocity stuff have not streamed
        pop = sim.populations[1, 2:2:end, 2:end]
        @fact pop - pop[1] => roughly(zeros(size(pop)...))
        pop = sim.populations[2, 1:2:end, 1:2:end]
        @fact pop - pop[1] => roughly(zeros(size(pop)...))
        @fact sim.populations[1, 1, 1] => not(roughly(sim.populations[1, 1, 2]))
    end

    context("no boundaries, just inlet") do
        sim.playground[:] = FLUID
        sim.playground[:, 1] = INLET
        sim.populations[:] = 0
        sim.populations[1, :, :] = 1
        sim.inlet_velocity = [1, 0]

        # Run one step
        run_lb(sim)
        fᵢ⁰ = sim.populations[:, 1, 2] 
        for i=1:size(sim.populations, 2), j=2:size(sim.populations, 2)
            @fact sim.populations[:, i, j] => roughly(fᵢ⁰)
        end
        for i=1:size(sim.populations, 2), j=2:size(sim.populations, 2)
            @fact sim.populations[:, i, 1] => not(roughly(fᵢ⁰))
        end
    end
	
	context("south no slip boundary, no inlet") do
        sim.playground[:] = FLUID
        sim.playground[end, :] = SOLID
        sim.populations[:] = 0
        sim.populations[5, :, end-1] = 1 # Downward velocity for particles next to wall

        # Run one step
        run_lb(sim)
        @fact sim.populations[3, :, end] => roughly(1)
    end
	
end

