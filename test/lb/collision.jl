using LatBo.LB: collision

facts("Collision kernels") do
    context("Single relaxation time") do
        context("Linear vs τ⁻¹") do
            @fact 2collision(1, [1], [0]) => collision(2, [1], [0])
            @fact 3collision(1, [0], [1]) => collision(3, [0], [1])
        end
        context("Linear vs fᵢ (and feq)") do
            f = rand(Int8, 10)
            @fact 2collision(1, f, zeros(f)) => collision(1, 2f, zeros(f))
            @fact 3collision(1, zeros(f), f) => collision(1, zeros(f), 3f)
        end
        context("Anti-symmetric vs fᵢ and feq") do
            fᵢ, feq = rand(Int8, 10), rand(Int8, 10)
            @fact collision(1, fᵢ, feq) => -collision(1, feq, fᵢ)
        end
        context("Lock-down the sign") do
            @fact collision(1, [1:10], zeros(Int8, 10)) .< 0 => all
        end
    end
end
