# Calculate the initial probability density function (we use the equation for the equilibrium function for this)

# Output will be all velocities for f1 at a particular lattice site -> f1[:,i,j]

using LatBo: equilibrium
using LatBo.D2Q9

function initial_probability(playground_var,initial_vel,initial_rho,e,weights,c) 

	if playground_var == 2
		return equilibrium([0.0; 0.0],e,weights,initial_rho,c)
	else
		return equilibrium(initial_vel,e,weights,initial_rho,c)
	end
end
