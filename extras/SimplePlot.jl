module SimplePlot
export plot
using LatBo.SandBox
using LatBo.Indices: gridcoords
using LatBo.LB: density, momentum
using LatBo.Playground
import Gadfly: plot
using Gadfly
function plot(what::Symbol, sim::SandBox, args...; kwargs...)
    x = [gridcoords(sim.indexing, i)[1] for i = 1:length(sim.playground)]
    y = [gridcoords(sim.indexing, i)[2] for i = 1:length(sim.playground)]
    const Npoints = size(sim.indexing)

    is_fluid = i -> sim.playground[i] == Playground.FLUID
    only_fluid = filter(is_fluid, collect(1:length(sim.playground)))
    nofluid = filter(i -> !is_fluid(i), collect(1:length(sim.playground)))

    const populations  = {parse("f_$i") => i for i in 1:9}

    title = Guide.title(string(what))
    if what == :playground
        plotthis = sim.playground
    elseif what == :density
        plotthis = density(sim)
    elseif what == :μx
        plotthis = momentum(sim)[1, :]
    elseif what == :μy
        plotthis = momentum(sim)[2, :]
    elseif what == :μ
        plotthis = sqrt(momentum(sim)[1, :].^2 + momentum(sim)[2, :].^2)
    elseif haskey(populations, what)
        j = populations[what]
        plotthis = reshape(sim.populations[j, :], Npoints)
        title = Guide.title("population: $(sim.lattice.celerities[:, j])")
    else
        error("Unknown quantity $what")
    end
    if what == :playground
        plot(x=x[only_fluid], y=y[only_fluid], color=plotthis[only_fluid],
            title, args...; kwargs...)
    else
        lines = reshape(plotthis[only_fluid], Npoints[1]-2, Npoints[2]-2)
        plot(
        layer(x=x[nofluid], y=y[nofluid], Theme(default_color=color("black")), Geom.point),
            layer(x=x[only_fluid], y=y[only_fluid], color=plotthis[only_fluid], Geom.point),
            layer(
                x=collect(2:(size(lines, 1)+1)), 
                y=collect(2:(size(lines, 2) + 1)),
                z=lines, Geom.contour
            ),
            title,
            args...; kwargs...
        )
    end
end
end
