using FactCheck: facts, context, @fact, not
using LatBo.geometry: is_in_pipe, is_in_half_space
using LatBo.Playground: SOLID, FLUID, INLET, OUTLET, NOTHING
import LatBo.Playground.initialize

facts("Constructing playground with a single pipe") do

    context("2d") do
        # Create playground a single pipe
        playpen = initialize((40, 40)) do i, j
            coords = Float64[i, j]
            if !is_in_pipe(coords, Float64[0, 1], Float64[20.5, 20.5], 18.)
                return SOLID
            elseif is_in_half_space(coords, Float64[0, -1], Float64[0, 1])
                return INLET
            elseif is_in_half_space(coords, Float64[0, 1], Float64[0, 40])
                return OUTLET
            else
                return FLUID
            end
        end

        # Check we got what we asked for
        @fact playpen .!= NOTHING --> all
        for indices in ([3, 20], [20, 5], [38, 35])
            @fact playpen[indices...] --> FLUID
        end
        for indices in ([2, 20], [39, 35], [2, 1], [2, 40])
            @fact playpen[indices...] --> SOLID
        end
        @fact playpen[3:end-3, 1] .== INLET --> all
        @fact playpen[3:end-3, 40] .== OUTLET --> all
    end

    context("3d") do
        # Create playground a single pipe
        playpen = initialize((10, 10, 20)) do i, j, k
            coords = Float64[i, j, k]
            if !is_in_pipe(coords, [0., 0, 1], [5.5, 5.5, 10.5], 4.)
                return SOLID
            elseif is_in_half_space(coords, [0., 0, -1], [0., 0, 1])
                return INLET
            elseif is_in_half_space(coords, [0., 0, 1], [0., 0, 20])
                return OUTLET
            else
                return FLUID
            end
        end

        # Check we got what we asked for
        @fact playpen .!= NOTHING --> all
        in_sphere(i, j) =  (i - 5.5)^2 + (j - 5.5)^2 <= 16
        for i = 1:size(playpen, 1), j = 1:size(playpen, 2)
            @fact playpen[i, j, 10] --> (in_sphere(i, j) ? FLUID: SOLID)
        end
        for i = 1:size(playpen, 1), j = 1:size(playpen, 2)
            @fact playpen[i, j, 1] --> (in_sphere(i, j) ? INLET: SOLID)
        end
        for i = 1:size(playpen, 1), j = 1:size(playpen, 2)
            @fact playpen[i, j, 20] --> (in_sphere(i, j) ? OUTLET: SOLID)
        end
    end
end
