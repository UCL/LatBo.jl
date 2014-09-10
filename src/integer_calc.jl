# This function calculates the new indices for assigning f2 in the
# streaming/propagation (can be easily extended to 3D)

function integer_calc(nx,ny,i,j,k,c)
	
	# Calculate the new integers to stream/propogate to
	newx = 1 + mod(i-1+c[1,k]+nx,nx)
	newy = 1 + mod(j-1+c[2,k]+ny,ny)
    return newx, newy
end
