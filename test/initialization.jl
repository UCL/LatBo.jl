using FactCheck: facts, @fact, roughly, @fact_throws
using LatBo: SingleRelaxationTime
using LatBo.playground: Feature, NOTHING

facts("Single relaxation time initialization") do
    data = SingleRelaxationTime(1.5, 1.6, 0.5, :D3Q15, (20, 20, 20))


    @fact data.Δt => roughly(1.5)
    @fact data.Δx => roughly(1.6)
    @fact data.τ⁻¹ => roughly(0.5)
    @fact data.kernel => :D3Q15

    @fact eltype(data.populations) => Float64
    @fact size(data.populations) => (20, 20, 20, 15)
    @fact data.populations => roughly(zeros(20, 20, 20, 15))

    @fact eltype(data.playground) => Feature
    @fact size(data.playground) => (20, 20, 20)
    @fact data.playground .== NOTHING => all
end

facts("Single relaxation time initialization failure modes") do
    @fact_throws SingleRelaxationTime(1.5, 1.6, 0.5, :D2Q15, (20, 20, 20))
    @fact_throws SingleRelaxationTime(1.5, 1.6, 0.5, :D4Q15, (20, 20, 20, 20))
    @fact_throws SingleRelaxationTime(1.5, 1.6, 0.5, :D1Q15, (20,))
    @fact_throws SingleRelaxationTime(-1.5, 1.6, 0.5, :D3Q15, (20, 20, 20))
    @fact_throws SingleRelaxationTime(1.5, -1.6, 0.5, :D3Q15, (20, 20, 20))
    @fact_throws SingleRelaxationTime(1.5, 1.6, -0.5, :D3Q15, (20, 20, 20))
end

