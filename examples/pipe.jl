using LatBo: SandBox, run!, Playground, geometry, Simulation
using LatBo.Indices: gridcoords
using LatBo.LB: FluidKernel, SingleRelaxationTime, ParabolicVelocityIOlet, NashZeroOrderPressure
using Gadfly
import LatBo.LB.velocity


println(1)
sim = SandBox(:D2Q9, (40, 40), ρ₀=1, ν₀=[1, 0])
println(2)
streamers = {
    Playground.INLET => ParabolicVelocityIOlet{Float64}([1, 0], [20, 20], 20, 1),
    Playground.OUTLET => ParabolicVelocityIOlet{Float64}([1, 0], [20, 20], 20, 0.4)
}

sim.kernels = {
    Playground.FLUID => FluidKernel(SingleRelaxationTime{Float64}(1.), streamers)
}
println(3)
Playground.initialize(sim.playground) do i, j
    coords = Float64[i, j]
    if !geometry.is_in_pipe(coords, Float64[0, 1], Float64[20.5, 20.5], 18.)
        return Playground.SOLID
    elseif geometry.is_in_half_space(coords, Float64[0, -1], Float64[0, 1])
        return Playground.INLET
    elseif geometry.is_in_half_space(coords, Float64[0, 1], Float64[0, 40])
        return Playground.OUTLET
    else
        return Playground.FLUID
    end
end


function velocity(sim::Simulation)
    velocities = zeros(size(sim.playground)..., 2)
    for (i, feature) in enumerate(sim.playground)
        if feature == Playground.FLUID
            velocities[gridcoords(sim.indexing, i)..., :] =
                velocity(sim.populations[:, i], sim.lattice.celerities)
        else
            velocities[gridcoords(sim.indexing, i)..., :] = [0, 0]
        end
    end
    velocities
end

println(4)
run!(sim, nsteps=20)
println(5)
alls = Any[]
push!(alls, velocity(sim))
npoints = size(sim.playground, 1)
plot(
    x=repeat([1:npoints...], inner=[npoints]),
    y=repeat([1:npoints...], outer=[npoints]),
    color=alls[1][:, :]
)
