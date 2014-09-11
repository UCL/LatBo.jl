using FactCheck: facts, context, @fact, not, roughly
using LatBo: noslip_boundary, integer_calc
using LatBo.playground: initialize
using LatBo.D2Q9

facts("No Slip Boundary Checks") do

	context("Box boundary") do
	
		grid = ones(Int32,4,3) # Initialize grid of FLUID points (x, y) == (i, j)
		# Specify solid points
			grid[:,1] = 2
			grid[:,size(grid,2)] = 2
			grid[1,:] = 2
			grid[size(grid,1),:] = 2
			
		f = zeros(Float64,9,4,3) # Initialize 2D lattice with all zero velocities
		f2 = zeros(Float64,9,4,3) # Initialize secondary grid
		c = D2Q9.celerities # Celerities Vector
		b_indices = [1 2 3 4 4 4 3 2 1 1;
					1 1 1 1 2 3 3 3 3 2] # Define boundary node indices as column vectors

		sol = transpose(Float64[0 0 0 0 0 8 0 8 0;
							0 0 5 0 5 17 0 17 0;
							0 0 14 0 14 0 9 0 9;
							0 0 0 0 0 0 18 0 18;
							0 11 0 11 0 0 0 0 0;
							0 0 0 0 0 15 0 15 0;
							0 0 12 0 12 6 0 6 0;
							0 0 3 0 3 0 16 0 16;
							0 0 0 0 0 0 7 0 7;
							0 4 0 4 0 0 0 0 0;])
		
		# All velocities of centre two points are unity
		fl_ind = [2 3; 2 2] # Fluid node indices as column vectors
		
		f[:,(fl_ind[:,1])...] = 1:9
		f[:,(fl_ind[:,2])...] = 10:18
		
		# Stream values to neighbour nodes
		for n = 1:2 # For each fluid node
			fl = fl_ind[:,n] # Read indices as column vector
			
			for k = 2:9 # For each velocity component find streamed indices
				s = integer_calc([size(grid)...],fl,c[:,k])
				
				# Generate new grid f2 with streamed values
				f2[k,(s)...] = f[k,(fl)...]
			end
		end
		f = f2 # Overwrite old values with new
		
		# For each boundary node apply boundary condition		
		for b = 1:size(b_indices,2)
			index = b_indices[:,b]
			noslip_boundary(grid,index,f)
			
			# Check the velocities match the solution
			println("Check boundary point ", b)
			@fact f[:,(index)...] => sol[:,b]
			
		end
		
    end
	
end
