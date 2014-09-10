using FactCheck: facts, context, @fact, not, roughly
using LatBo

facts("Check that streaming integers are correctly calculated") do

# Can freely edit these inputs (although i and j must not be greater than nx or ny)
nx = 4
ny = 3
i = 4
j = 3

# These inputs must not be changed
nVels = 9
c = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1]
f1 = zeros(nVels,nx,ny)
f2 = zeros(nVels,nx,ny)

# Assign 1.0 to all 9 velocities at (i,j) lattice site
for k = 1:nVels
	f1[k,i,j] = 1.0
end

# Call integer_calc and assign new f2 values from f1
for k = 1:nVels
	
	# Calculate the appropriate integers for f2
	newx = LatBo.integer_calc(nx,i,c[1,k])
	newy = LatBo.integer_calc(ny,j,c[2,k])

	# Assign 
	f2[k,newx,newy] = f1[k,i,j]
end

# Set new variables for testing
i_plus = i + 1
i_minus = i - 1
j_plus = j + 1
j_minus = j - 1

# Conditional statements in case i or j are on a boundary
if i_plus == nx + 1
	i_plus = 1
end
if i_minus == 0
	i_minus = nx
end
if j_plus == ny + 1
	j_plus = 1
end
if j_minus == 0
	j_minus = ny
end

# Check the velocities are moving in the right direction (ie. the correct integers are calculated by integer_calc)
@fact f2[1,i,j] => roughly(1.0)
@fact f2[2,i_plus,j] => roughly(1.0)
@fact f2[3,i,j_plus] => roughly(1.0)
@fact f2[4,i_minus,j] => roughly(1.0)
@fact f2[5,i,j_minus] => roughly(1.0)
@fact f2[6,i_plus,j_plus] => roughly(1.0)
@fact f2[7,i_minus,j_plus] => roughly(1.0)
@fact f2[8,i_minus,j_minus] => roughly(1.0)
@fact f2[9,i_plus,j_minus] => roughly(1.0)

end
