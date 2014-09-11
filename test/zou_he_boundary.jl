using FactCheck: facts, context, @fact, not, roughly
using LatBo

facts("Check that Zou-He Boundaries are correctly calculated") do

# Can freely edit these inputs
nx = 5
ny = 7

# Initialise boundary velocity (NOTE: second integer must be consistent with the boundary you are calculating)
bound_vel = zeros(2,nx)

# Assign this array to the unknown lattice velocities (eg 2,6,9 for west boundary)
#unknown_vels = [2 6 9]	# West
#unknown_vels = [4 8 7]	# East
#unknown_vels = [3 6 7]	# South
unknown_vels = [5 8 9] # North

# Set the boundary velocity array
for index = 2:nx-1
	bound_vel[1,index] = 2.0
	bound_vel[2,index] = 2.5
end

# Assign what values our test variables will take
c1 = 3.2
c2 = 1/9
c3 = 10/3
c4 = 15.0
c5 = 19.6667
c6 = 11.2223
c7 = 12000
c8 = 7.0
c9 = 5.5

# These inputs must remain the same
nVels = 9
f2 = zeros(nVels,nx,ny)

# Assign test values to f2
for i = 2:nx-1
	for j = ny:ny

		# Assign index to the loop variable (saves me having to change it everytime I want to look at N,S or E,W boundaries
		index = i

		f2[1,i,j] = c1
		f2[2,i,j] = c2
		f2[3,i,j] = c3
		f2[4,i,j] = c4
		f2[5,i,j] = c5
		f2[6,i,j] = c6
		f2[7,i,j] = c7
		f2[8,i,j] = c8
		f2[9,i,j] = c9

		# Reset the unknown values to 0.0
		f2[unknown_vels,i,j] = 0.0


		# Call the function to be tested
		f2[:,i,j] = LatBo.zou_he_boundary(i,j,nx,ny,f2[:,i,j],bound_vel[:,index])

		# Test the correct values are calculated
		# West boundary
		if i == 1

			@fact f2[2,i,j] => roughly(c4 + bound_vel[1,index] * 2.0 / 3.0)
			@fact f2[6,i,j] => roughly(c8 + bound_vel[1,index] * 1.0 / 6.0 + 0.5 * (c5 - c3) + 0.5 * bound_vel[2,index])
			@fact f2[9,i,j] => roughly(c7 + bound_vel[1,index] * 1.0 / 6.0 - 0.5 * (c5 - c3) - 0.5 * bound_vel[2,index])

		# East boundary
		elseif i == nx

			@fact f2[4,i,j] => roughly(c2 - bound_vel[1,index] * 2.0 / 3.0)
			@fact f2[8,i,j] => roughly(c6 - bound_vel[1,index] * 1.0 / 6.0 - 0.5 * (c5 - c3) - 0.5 * bound_vel[2,index])
			@fact f2[7,i,j] => roughly(c9 - bound_vel[1,index] * 1.0 / 6.0 + 0.5 * (c5 - c3) + 0.5 * bound_vel[2,index])

		# South boundary
		elseif j == 1

			@fact f2[3,i,j] => roughly(c5 + bound_vel[2,index] * 2.0 / 3.0)
			@fact f2[6,i,j] => roughly(c8 + bound_vel[2,index] * 1.0 / 6.0 + 0.5 * (c4 - c2) + 0.5 * bound_vel[1,index])
			@fact f2[7,i,j] => roughly(c9 + bound_vel[2,index] * 1.0 / 6.0 - 0.5 * (c4 - c2) - 0.5 * bound_vel[1,index])

		# North boundary
		elseif j == ny
			
			@fact f2[5,i,j] => roughly(c3 - bound_vel[2,index] * 2.0 / 3.0)
			@fact f2[8,i,j] => roughly(c6 - bound_vel[2,index] * 1.0 / 6.0 - 0.5 * (c4 - c2) - 0.5 * bound_vel[1,index])
			@fact f2[9,i,j] => roughly(c7 - bound_vel[2,index] * 1.0 / 6.0 + 0.5 * (c4 - c2) + 0.5 * bound_vel[1,index])
		end
	end
end
end
