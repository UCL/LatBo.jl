function equilibrium(u,e,weights,rho,c)
#= Calculation of incompressible collision operator for a given lattice point

u 			= Macroscopic velocity vector
e			= Unit velocity matrix (celerities)
weights		= Weights of Gauss-Hermite quadrature
rho 			= Density
c			= Sound speed

Equilibrium probability function for i-th direction feq_i of the form:

weight_i * [rho + A_i + B_i - C]

=#

# Declare empty feq array with length based on number of velocities found from length of weights
feq = zeros(Float64,1,size(weights,2))

	# Calculate C
	C = dot(u,u) / (2*c^2)

	# Get number of velocities from size of weights and loop through
	for i = 1:size(weights,2)
		
		# Calculate A
		A = dot(e[:,i],u) / c^2
		
		# Calculate B
		B = (dot(e[:,i],u))^2 / (2*c^4)
		
		# Calculate feq_i
		feq[i] = weights[i] * (rho + A + B - C)
		
	end
	
	return feq
		
end # function equilibrium
 

############################################################
function collision(u,f,e,weights,rho,τ⁻¹,c)
#= Calculate collision operator at a given lattice point using equilibrium distribution feq computed using equilibrium()

u 			= Macroscopic velocity vector
f 			= Probability distribution function
e			= Unit velocity matrix (celerities)
weights		= Weights of Gauss-Hermite quadrature
rho 			= Density
τ⁻¹ 			= Relaxation frequency
c			= Sound speed

=#

τ⁻¹ * (equilibrium(u,e,weights,rho,c) - f)

end # function collision
