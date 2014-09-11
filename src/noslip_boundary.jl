using LatBo: D2Q9
using LatBo.playground: SOLID, FLUID

function noslip_boundary(grid, index, f, mid=false)

# grid		= array indicating type of point (FLUID, SOLID , ... etc)
# index		= vector of indices corresponding to particular boundary point
# f			= probability distribution function for grid
# mid		= set TRUE to use midway boundary condition, else will use full bounce back.

lattice_size	= [(size(grid))...]

directions = D2Q9.celerities
m_opp = Int32[0 3 4 1 2 7 8 5 6] # Opposite directions

# Consider the lattice point (j,k)
for m = 2:size(directions,2) # Only consider the neighbours

	direction = directions[:,m] # Consider each direction

	if any( (index + direction) .> lattice_size)
		# Outside the lattice
		# println("Outside +Boundary")
	
	elseif any( (index + direction) .< ones(size(index)...) )
		# Also outside lattice
		# println("Outside -Boundary")
	
	elseif grid[(index + direction)...] == FLUID
		if mid
			f[m,(index)...] = f[m_opp[m]+1,(index+direction)...]
		else
			f[m,(index)...] = f[m_opp[m]+1,(index)...]
		end
	end
end
	

end