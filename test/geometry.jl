using FactCheck: facts, context, @fact
using LatBo.geometry: is_inside_pipe, is_in_half_space, is_in_sphere

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
            @fact is_inside_pipe(r, d₀, center, radius) => (y <= radius)
        end
    end

    context("3d") do
        # Direction
        d₀ = Float64[1, 1, 1]
        d⟂⁰ = Float64[0, -1, 1]/√2
        d⟂¹ = Float64[-2, 1, 1]/√(4+1+1)
        center = Float64[5, 0, 0]
        radius = 3.

        yₛ = radius * [rand(N)..., (1. + rand(N))...]
        zₛ = radius * [rand(N)..., (1. + rand(N))...]
        for x in rand(N) * 4., y in yₛ, z in zₛ
            r = center + x * d₀ + y * d⟂⁰ + z * d⟂¹
            inside = is_inside_pipe
            @fact inside(r, d₀, center, radius) => (norm([z, y]) <= radius)
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

        xₛ = 4.0*[rand(N)..., -rand(N)...]
        yₛ = 4.0(rand(N) - 0.5)
        zₛ = 4.0(rand(N) - 0.5)
        for x in xₛ, y in yₛ, z in zₛ
            r = Ω₀ + x * d₀ + y * d⟂⁰ + z * d⟂¹
            inside = is_in_half_space
            @fact inside(r, d₀, Ω₀) => (x >= 0)
        end
    end
end

facts("Check sphere geometries") do
    N = 4  # Measure of the number of points to test
    context("3d") do
        # Direction
        Ω₀ = Float64[4, 3, 1]
        radius = 3.

        xₛ = 0.5radius*[rand(N)..., 1. + rand(N)...]
        yₛ = 0.5radius*[rand(N)..., 1. + rand(N)...]
        zₛ = 0.5radius*[rand(N)..., 1. + rand(N)...]
        for x in xₛ, y in yₛ, z in zₛ
            r = [x, y, z] + Ω₀
            inside = is_in_sphere
            @fact inside(r, Ω₀, radius) => (norm([x, y, z]) <= radius)
        end
    end
end
