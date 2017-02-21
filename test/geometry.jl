using LatBo.geometry: is_in_pipe, is_in_half_space, is_in_sphere

facts("Check pipe geometries") do
    N = 4  # Measure of the number of points to test

    context("2d") do
        # Direction
        d₀ = Float64[1, 1]
        d⟂ = Float64[1, -1]/√2
        center = Float64[5, 0]
        radius = 3.

        for x in rand(N) * 4., y in radius * [rand(N)..., (1. + rand(N))...]
            r = center + x * d₀ + y * d⟂
            @fact is_in_pipe(r, d₀, center, radius) --> (y <= radius)
        end
    end

    context("3d") do
        # Direction
        d₀ = Float64[1, 1, 1]
        d⟂⁰ = Float64[0, -1, 1]/√2
        d⟂¹ = Float64[-2, 1, 1]/√(4+1+1)
        center = Float64[5, 0, 0]
        radius = 3.
        in_pipe(r) = is_in_pipe(r, d₀, center, radius)

        ys = radius * [rand(N)...]
        zs = radius * [rand(N)..., (1. + rand(N))...]
        for x in rand(N) * 4., y in ys, z in zs
            if norm([z, y]) <= radius
                @fact center + x * d₀ + y * d⟂⁰ + z * d⟂¹ --> in_pipe
            else
                @fact center + x * d₀ + y * d⟂⁰ + z * d⟂¹ --> not(in_pipe)
            end
        end
    end
end

facts("Check half-plane geometries") do
    N = 4  # Measure of the number of points to test
    context("3d") do
        # Direction
        d₀ = Float64[-1, -1, -1]
        d⟂⁰ = Float64[0, -1, 1]/√2
        d⟂¹ = Float64[-2, 1, 1]/√(4+1+1)
        Ω₀ = Float64[4, 3, 1]

        xs = 4.0*[rand(N)..., -rand(N)...]
        ys = 4.0(rand(N) - 0.5)
        zs = 4.0(rand(N) - 0.5)
        for x in xs, y in ys, z in zs
            r = Ω₀ + x * d₀ + y * d⟂⁰ + z * d⟂¹
            inside = is_in_half_space
            @fact inside(r, d₀, Ω₀) --> (x >= 0)
        end
    end
end

facts("Check sphere geometries") do
    N = 4  # Measure of the number of points to test
    context("3d") do
        # Direction
        Ω₀ = Float64[4, 3, 1]
        radius = 3.

        xs = 0.5radius*[rand(N)..., 1. + rand(N)...]
        ys = 0.5radius*[rand(N)..., 1. + rand(N)...]
        zs = 0.5radius*[rand(N)..., 1. + rand(N)...]
        for x in xs, y in ys, z in zs
            r = [x, y, z] + Ω₀
            inside = is_in_sphere
            @fact inside(r, Ω₀, radius) --> (norm([x, y, z]) <= radius)
        end
    end
end
