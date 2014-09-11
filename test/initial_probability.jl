using FactCheck: facts, context, @fact, not, roughly
using LatBo

facts("Check that inital_probability works correctly") do

# Set the initial velocity and density
initial_vel = [15.0; -1/3]
initial_rho = 1.5

# Constant parameters
nVels = 9
weights = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]
e = D2Q9.celerities
c = 1.0 / sqrt(3.0)

# Flags for feature type at each lattice site
nx = 2
ny = 2
playground_var = [1 2;
				3 4]

# Initialise the probability distribution function
f1 = zeros(nVels,nx,ny)

# Loop through each lattice site and call initial_probability
for i = 1:nx
	for j = 1:ny

		f1[:,i,j] = LatBo.initial_probability(playground_var[i,j],initial_vel,initial_rho,e,weights,c)

	end
end


@fact initial_rho => roughly(sum(f1[:,1,1]))
@fact initial_vel => roughly(e * f1[:,1,1])

@fact initial_rho => roughly(sum(f1[:,1,2]))
@fact [0.0; 0.0] => roughly(e * f1[:,1,2])

@fact initial_rho => roughly(sum(f1[:,2,1]))
@fact initial_vel => roughly(e * f1[:,2,1])

@fact initial_rho => roughly(sum(f1[:,2,2]))
@fact initial_vel => roughly(e * f1[:,2,2])

end

