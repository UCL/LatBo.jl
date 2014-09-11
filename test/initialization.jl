using FactCheck: facts, @fact, roughly, @fact_throws
using LatBo: SingleRelaxationTime
using LatBo.playground: Feature, NOTHING

facts("Single relaxation time initialization") do
    data = SingleRelaxationTime(0.5, (20, 20, 20))


    @fact data.τ⁻¹ => roughly(0.5)

    @fact eltype(data.populations) => Float64
    @fact size(data.populations) => (19, 20, 20, 20)
    @fact data.populations => roughly(zeros(19, 20, 20, 20))

    @fact eltype(data.playground) => Feature
    @fact size(data.playground) => (20, 20, 20)
    @fact data.playground .== NOTHING => all
end

facts("Single relaxation time initialization failure modes") do
    @fact_throws SingleRelaxationTime(0.5, (20, 20, 20, 20))
    @fact_throws SingleRelaxationTime(0.5, (20,))
end

