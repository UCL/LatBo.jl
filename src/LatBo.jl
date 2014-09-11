module LatBo

export geometry, playground, LatticeBoltzmann, SingleRelaxationTime, D2Q9,
    D3Q19, thermodynamics, collision, lattice_loop, integer_calc,
	noslip_boundary, run_lb#, visualisation

abstract LatticeBoltzmann

include("geometry.jl")
include("playground.jl")
include("single_relaxation_time.jl")
#include("plot_frame.jl")
include("thermodynamics.jl")
include("collision.jl")
include("integer_calc.jl")
include("kernel.jl")
include("zou_he_boundary.jl")
include("initial_probability.jl")
include("noslip_boundary.jl")
#include("visualisation.jl")


# Runs lattice boltzmann for n steps
function run_lb(observer::Function, sim::LatticeBoltzmann, nsteps::Int)
    for n = 1:nsteps
        run_lb(observer, sim)
    end
end

# Runs lattice boltzmann for single step
function run_lb(observer::Function, sim::LatticeBoltzmann)
    # Aliases for easier acces to quantities
    this_pop = sim.populations
    next_pop = sim.next_populations
    gridsize = [size(sim.playground)...]
    kernel = sim.kernel
    celerities = sim.kernel.celerities

    # Loop over each lattice site
    lattice_loop(sim) do indices, fᵢ, feature
        # Apply collision step
		this_pop[:, indices...] += collision(fᵢ, kernel, sim.τ⁻¹)
        # Apply streaming step
        for v = 1:size(celerities, 2)
            streamed = integer_calc(gridsize, indices, celerities[:, v])
            next_pop[v, streamed...] = this_pop[v, indices...]
        end
    end

    # Applying boundary condition
    lattice_loop(sim) do indices, fᵢ, feature
        if feature == playground.INLET
            next_pop[:, indices...] = zou_he_boundary(
                indices..., gridsize..., fᵢ, sim.inlet_velocity)
        elseif feature == playground.SOLID
            noslip_boundary(sim.playground,indices,next_pop)
        end
    end

    # reassign populations and next_populations with swap
    sim.populations  = next_pop
    sim.next_populations = this_pop

    # run observer at each step
    observer()
end

# Run simulation for N numbers of steps
run_lb(sim::LatticeBoltzmann, nsteps::Int) = (
    run_lb(()->nothing, sim, nsteps)
)
run_lb(sim::LatticeBoltzmann) = run_lb(()->nothing, sim)

end # module
