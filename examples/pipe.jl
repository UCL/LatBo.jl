using LatBo
using Gadfly


sim = SingleRelaxationTime(1., (20, 20))
sim.inlet_velocity = [.1, 0]
playground.initialize(sim.playground) do i, j
    coords = Float64[i, j]
    if !geometry.is_in_pipe(coords, Float64[0, 1], Float64[10.5, 10.5], 8.)
        return playground.SOLID
    elseif geometry.is_in_half_space(coords, Float64[0, -1], Float64[0, 1])
        return playground.INLET
    elseif geometry.is_in_half_space(coords, Float64[0, 1], Float64[0, 20])
        return playground.OUTLET
    else
        return playground.FLUID
    end
end

lattice_loop(sim) do indices, populations, feature
    sim.populations[:, indices...] = LatBo.initial_probability(
       feature, [.1, 0], 5.,
       sim.kernel.celerities, sim.kernel.weights, sim.kernel.speed_of_sound
    )
end

function velocity(sim::LatticeBoltzmann)
    velocities = zeros(size(sim.playground)..., 2)
    lattice_loop(sim) do indices, fᵢ, feature
        if feature != playground.SOLID
            velocities[indices..., :] = thermodynamics.velocity(
                fᵢ, sim.kernel.D2Q9.celerities)
        else
            velocities[indices..., :] = [-0, -0]
        end
    end
    velocities
end

alls = Any[]
push!(alls, velocity(sim))
LatBo.run_lb(sim, 5) do
    push!(alls, velocity(sim))
end
npoints = size(sim.playground, 1)
plot(
    x=repeat([1:npoints...], inner=[npoints]),
    y=repeat([1:npoints...], outer=[npoints]),
    color=alls[end]
)
