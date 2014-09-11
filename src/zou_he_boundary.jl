# Apply Zou-He boundary conditions

function zou_he_boundary(i,j,nx,ny,f2,bound_vel)

	# Determine where the Boundary is (N,S,E,W)
	# Boundary is on the west boundary	
	if i == 1

		f2[2] = f2[4] + bound_vel[1] * 2.0 / 3.0
		f2[6] = f2[8] + bound_vel[1] * 1.0 / 6.0 + 0.5 * (f2[5] - f2[3]) + 0.5 * bound_vel[2]
		f2[9] = f2[7] + bound_vel[1] * 1.0 / 6.0 - 0.5 * (f2[5] - f2[3]) - 0.5 * bound_vel[2]
	
	# Boundary is on the east boundary
	elseif i == nx

		f2[4] = f2[2] - bound_vel[1] * 2.0 / 3.0
		f2[8] = f2[6] - bound_vel[1] * 1.0 / 6.0 - 0.5 * (f2[5] - f2[3]) - 0.5 * bound_vel[2]	
		f2[7] = f2[9] - bound_vel[1] * 1.0 / 6.0 + 0.5 * (f2[5] - f2[3]) + 0.5 * bound_vel[2]

	
	# Boundary is on the south boundary
	elseif j == 1

		f2[3] = f2[5] + bound_vel[2] * 2.0 / 3.0
		f2[6] = f2[8] + bound_vel[2] * 1.0 / 6.0 + 0.5 * (f2[4] - f2[2]) + 0.5 * bound_vel[1]
		f2[7] = f2[9] + bound_vel[2] * 1.0 / 6.0 - 0.5 * (f2[4] - f2[2]) - 0.5 * bound_vel[1]
	
	# Boundary is on the north boundary
	elseif j == ny

		f2[5] = f2[3] - bound_vel[2] * 2.0 / 3.0		
		f2[8] = f2[6] - bound_vel[2] * 1.0 / 6.0 - 0.5 * (f2[4] - f2[2]) - 0.5 * bound_vel[1]
		f2[9] = f2[7] - bound_vel[2] * 1.0 / 6.0 + 0.5 * (f2[4] - f2[2]) + 0.5 * bound_vel[1]
	end

	return f2

end
