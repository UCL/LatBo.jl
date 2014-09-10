# Incompressible channel flow with Zou-He boundary conditions on inlet, outlet and walls
# Joe O'Connor
# 09/09/14

#### FUNCTION DEFINITIONS ####
## Function for the collision and stream steps
function collide_stream(nx,ny,nVels,u,v,u_eq,v_eq,rho,cx,cy,w,feq,f1,f2)
	
	# Loop through all nodes
	for i = 1:nx
		for j = 1:ny
			
			# Assign equilibrium velocity (force terms can be included)
			u_eq[i,j] = u[i,j]
			v_eq[i,j] = v[i,j]

			# Calculate first coefficient for feq
			C1 = u_eq[i,j] * u_eq[i,j] + v_eq[i,j] * v_eq[i,j]

			# Loop through all velocities
			for k = 1:nVels
				
				# Calculate second coefficient for feq
				C2 = u_eq[i,j] * cx[k] + v_eq[i,j] * cy[k]

				# Calculate equilibrium function
				feq[k] = w[k] * (rho[i,j] + 3.0 * C2 + 4.5 * C2 * C2 - 1.5 * C1)

				# Collision step
				f1[k,i,j] = f1[k,i,j] * (1.0 - omega) + feq[k] * omega

				# Streaming step
				# Check if integer
				newx = 1 + mod(i-1+cx[k]+nx,nx)
				newy = 1 + mod(j-1+cy[k]+ny,ny)
				f2[k,newx,newy] = f1[k,i,j]
			end
		end
	end
	return f2
end

## Function for setting the boundary conditions
function boundaries(nx,ny,u,v,rho,f2)

	# Zou-He velocity BCs on all boundaries
	# Top and bottom (no-slip) walls
	for i = 2:nx-1
		
		# Bottom wall
		j = 1
		u[i,j] = 0.0
		v[i,j] = 0.0
		rho[i,j] = v[i,j] + f2[1,i,j] + f2[2,i,j] + f2[4,i,j] + 2.0 * (f2[5,i,j] + f2[8,i,j] + f2[9,i,j])

		f2[3,i,j] = f2[5,i,j] + v[i,j] * 2.0 / 3.0
		f2[6,i,j] = f2[8,i,j] + v[i,j] * 1.0 / 6.0 + 0.5 * (f2[4,i,j] - f2[2,i,j]) + 0.5 * u[i,j]
		f2[7,i,j] = f2[9,i,j] + v[i,j] * 1.0 / 6.0 - 0.5 * (f2[4,i,j] - f2[2,i,j]) - 0.5 * u[i,j]

		# Top wall
		j = ny
		u[i,j] = 0.0
		v[i,j] = 0.0
		rho[i,j] = -v[i,j] + f2[1,i,j] + f2[2,i,j] + f2[4,i,j] + 2.0 * (f2[3,i,j] + f2[6,i,j] + f2[7,i,j])

		f2[5,i,j] = f2[3,i,j] - v[i,j] * 2.0 / 3.0
		f2[8,i,j] = f2[6,i,j] - v[i,j] * 1.0 / 6.0 + 0.5 * (f2[2,i,j] - f2[4,i,j]) - 0.5 * u[i,j]
		f2[9,i,j] = f2[7,i,j] - v[i,j] * 1.0 / 6.0 - 0.5 * (f2[2,i,j] - f2[4,i,j]) + 0.5 * u[i,j]
	end

	# Inlet and outlet
	for j = 2:ny-1

		# Inlet
		i = 1
		u[i,j] = u_avg
		v[i,j] = 0.0
		rho[i,j] = u[i,j] + f2[1,i,j] + f2[3,i,j] + f2[5,i,j] + 2.0 * (f2[4,i,j] + f2[7,i,j] + f2[8,i,j])

		f2[2,i,j] = f2[4,i,j] + u[i,j] * 2.0 / 3.0
		f2[6,i,j] = f2[8,i,j] + u[i,j] * 1.0 / 6.0 + 0.5 * (f2[5,i,j] - f2[3,i,j]) + 0.5 * v[i,j]
		f2[9,i,j] = f2[7,i,j] + u[i,j] * 1.0 / 6.0 - 0.5 * (f2[5,i,j] - f2[3,i,j]) - 0.5 * v[i,j]

		# Outlet
		i = nx
		u[i,j] = u[i-1,j]
		v[i,j] = 0.0
		rho[i,j] = -u[i,j] + f2[1,i,j] + f2[3,i,j] + f2[5,i,j] + 2.0 * (f2[2,i,j] + f2[6,i,j] + f2[9,i,j])

		f2[4,i,j] = f2[2,i,j] - u[i,j] * 2.0 / 3.0
		f2[7,i,j] = f2[9,i,j] - u[i,j] * 1.0 / 6.0 + 0.5 * (f2[5,i,j] - f2[3,i,j]) + 0.5 * v[i,j]
		f2[8,i,j] = f2[6,i,j] - u[i,j] * 1.0 / 6.0 - 0.5 * (f2[5,i,j] - f2[3,i,j]) - 0.5 * v[i,j]
	end

	# Corner boundaries
	# Bottom left corner
	i = 1
	j = 1
	rho[i,j] = rho[i,j+1]
	u[i,j] = 0.0
	v[i,j] = 0.0
	f2[2,i,j] = f2[4,i,j] + u[i,j] * 2.0 / 3.0
	f2[3,i,j] = f2[5,i,j] + v[i,j] * 2.0 / 3.0
	f2[6,i,j] = f2[8,i,j] + (u[i,j] + v[i,j]) * 1.0 / 6.0
	f2[7,i,j] = (-u[i,j] + v[i,j]) * 1.0 / 12.0
	f2[9,i,j] = (u[i,j] - v[i,j]) * 1.0 / 12.0
	f2[1,i,j] = 0.0
	vsum = f2[1:end,i,j]
	vsum = sum(vsum)
	f2[1,i,j] = rho[i,j] - vsum

	# Top left corner
	i = 1
	j = ny
	rho[i,j] = rho[i,j-1]
	u[i,j] = 0.0
	v[i,j] = 0.0
	f2[2,i,j] = f2[4,i,j] + u[i,j] * 2.0 / 3.0
	f2[5,i,j] = f2[3,i,j] - v[i,j] * 2.0 / 3.0
	f2[9,i,j] = f2[7,i,j] + (u[i,j] - v[i,j]) * 1.0 / 6.0
	f2[6,i,j] = (u[i,j] + v[i,j]) * 1.0 / 12.0
	f2[8,i,j] = (-u[i,j] - v[i,j]) * 1.0 / 12.0
	f2[1,i,j] = 0.0
	vsum = f2[1:end,i,j]
	vsum = sum(vsum)
	f2[1,i,j] = rho[i,j] - vsum

	# Bottom right corner
	i = nx
	j = 1
	rho[i,j] = rho[i,j+1]
	u[i,j] = 0.0
	v[i,j] = 0.0
	f2[4,i,j] = f2[2,i,j] - u[i,j] * 2.0 / 3.0
	f2[3,i,j] = f2[5,i,j] + v[i,j] * 2.0 / 3.0
	f2[7,i,j] = f2[9,i,j] + (-u[i,j] + v[i,j]) * 1.0 / 6.0
	f2[6,i,j] = (u[i,j] + v[i,j]) * 1.0 / 12.0
	f2[8,i,j] = (-u[i,j] - v[i,j]) * 1.0 / 12.0
	f2[1,i,j] = 0.0
	vsum = f2[1:end,i,j]
	vsum = sum(vsum)
	f2[1,i,j] = rho[i,j] - vsum

	# Top right corner
	i = nx
	j = ny
	rho[i,j] = rho[i,j-1]
	u[i,j] = 0.0
	v[i,j] = 0.0
	f2[4,i,j] = f2[2,i,j] - u[i,j] * 2.0 / 3.0
	f2[5,i,j] = f2[3,i,j] - v[i,j] * 2.0 / 3.0
	f2[8,i,j] = f2[6,i,j] + (-u[i,j] - v[i,j]) * 1.0 / 6.0
	f2[7,i,j] = (-u[i,j] + v[i,j]) * 1.0 / 12.0
	f2[9,i,j] = (u[i,j] - v[i,j]) * 1.0 / 12.0
	f2[1,i,j] = 0.0
	vsum = f2[1:end,i,j]
	vsum = sum(vsum)
	f2[1,i,j] = rho[i,j] - vsum

	# Return the appropriate variables to caller
	return u, v, rho, f2
end

## Function for calculating rho, u and v
function rho_u_v(nx,ny,nVels,cx,cy,rho,u,v,f2)

	# Calculate rho, u and v
	for i = 1:nx
		for j = 1:ny

			# Reset to zero for summing the moments
			rho[i,j] = 0.0
			u[i,j] = 0.0
			v[i,j] = 0.0
	
			# Sum through all velocities
			for k = 1:nVels
				rho[i,j] = rho[i,j] + f2[k,i,j]
				u[i,j] = u[i,j] + cx[k] * f2[k,i,j]
				v[i,j] = v[i,j] + cy[k] * f2[k,i,j]
			end
		end
	end
	return rho, u, v
end

## Function for writing the values to data files
function data_writer(u,v,rho)
	
	# Write values to data files for MATLAB plot
	f = open("DATA/u_vel.dat","w")
	write(f,"$u")
	close(f)
	
	f = open("DATA/v_vel.dat","w")
	write(f,"$v")
	close(f)
	
	f = open("DATA/rho.dat","w")
	write(f,"$rho")
	close(f)
	
	f = open("DATA/x_y.dat","w")
	write(f,"$x_long\n$y_long")
	close(f)
end



#### BEGINNING OF MAIN CODE ####
## SETUP THE PROBLEM PARAMETERS
# Lattice parameters (D2Q9)
tic()
nVels = 9
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]
cx = [0 1 0 -1 0 1 -1 -1 1]
cy = [0 0 1 0 -1 1 1 -1 -1]


# Problem setup
Re = 12.0
u_max = 0.2
tau = 1.0
aspect_ratio = 5.0
develop_coeff = 2.0

# Problem calculations
u_avg = u_max * 2.0 / 3.0
omega = 1.0 / tau
k_visc = (tau - 0.5) / 3.0
y_long = iround(Re * k_visc / u_avg)
x_long = iround(y_long * aspect_ratio)
nx = x_long + 1
ny = y_long + 1
nSteps = iround((x_long / u_avg) * develop_coeff)

# Print domain size
@printf "Domain size is %i X %i\n" nx ny

# Initial (uniform) conditions for rho, u and v
rho_0 = 5.0
u_0 = u_avg
v_0 = 0.0

# Initialise the macroscopic matrices
rho = zeros(Float64,nx,ny)
u = zeros(Float64,nx,ny)
v = zeros(Float64,nx,ny)

# Fill the matrices with the initial conditions
for i = 1:nx
	
	# Initial density field
	for j = 1:ny
		rho[i,j] = rho_0
	end
	
	# Initial velocity field (u=v=0 at the wall)
	for j = 2:ny-1
		u[i,j] = u_0
		v[i,j] = v_0
	end
end

# Initalise distribution functions and equilibrium velocities (for force terms) to zero
u_eq = zeros(Float64,nx,ny)
v_eq = zeros(Float64,nx,ny)
feq = zeros(Float64,nVels)
f1 = zeros(Float64,nVels,nx,ny)
f2 = zeros(Float64,nVels,nx,ny)



## CALCULATE F1 AND FEQ BASED ON INITIAL CONDITIONS
for i = 1:nx
	for j = 1:ny

		# Assign equilibrium velocities (in case we want to include extra force terms later)
		u_eq[i,j] = u[i,j]
		v_eq[i,j] = v[i,j]

		# Calculate first coefficient for feq
		# Might need to get rid of the -0
		C1 = u_eq[i,j] * u_eq[i,j] + v_eq[i,j] * v_eq[i,j]

		# Loop through all velocities
		for k = 1:nVels
			
			# Calculate second coefficient for feq
			# Might need to get rid of the -0
			C2 = u_eq[i,j] * cx[k] + v_eq[i,j] * cy[k]

			# Calculate equilibrium function for feq
			feq[k] = w[k] * (rho[i,j] + 3.0 * C2 + 4.5 * C2 * C2 - 1.5 * C1)

			# Initially set f1 to feq
			f1[k,i,j] = feq[k]
		end
	end
end


## MAIN LBM ALGORITHM
for counter = 1:nSteps

	# Collide and stream
	f2 = collide_stream(nx,ny,nVels,u,v,u_eq,v_eq,rho,cx,cy,w,feq,f1,f2)

	# Set boundary conditions
	u, v, rho, f2 = boundaries(nx,ny,u,v,rho,f2)

	# Calculate new rho, u and v
	rho, u, v = rho_u_v(nx,ny,nVels,cx,cy,rho,u,v,f2)

	# Reassign f1 to f2
	for i = 1:nx
		for j = 1:ny
			for k = 1:nVels
				f1[k,i,j] = f2[k,i,j]
			end
		end
	end

	# Print timestep information
	@printf "Finished timestep %i of %i\n" counter nSteps

end

data_writer(u,v,rho)
toc()
