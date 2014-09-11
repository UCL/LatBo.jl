using FactCheck: facts, context, @fact, not, roughly
using LatBo: collision, equilibrium

facts("Check collision function") do

    context("Zero test") do
		u = zeros(Float64,2)
		f = zeros(Float64,10)
		e = zeros(Float64,2,10)
		weights = zeros(Float64,10)
		rho = 0
		τ⁻¹ = 0
		c = 1
       
        @fact collision(u,f,e,weights,rho,τ⁻¹,c) => roughly(zeros(Float64, length(weights)))
    end
	
	context("Ones test") do
		u = ones(Float64,2)
		f = ones(Float64,10)
		e = ones(Float64,2,10)
		weights = ones(Float64,10)
		rho = 1
		τ⁻¹ = 1
		c = 1
       
        @fact collision(u,f,e,weights,rho,τ⁻¹,c) => roughly(3*ones(Float64,length(weights)))
    end

end
