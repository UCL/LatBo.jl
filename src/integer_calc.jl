# This function calculates the new indices for assigning f2 in the
# streaming/propagation (can be easily extended to 3D)

function integer_calc(nLattice,index,velocity)
	
	# Calculate the new integers to stream/propogate to
	1 + mod(index-1+velocity,nLattice)  
end
