function  open_outlet(indices,gridsize,f2)

	# If outlet is on west boundary
	if indices[1] == 1
		return f2[:,indices[1]+1,indices[2]]

	# If outlet is on east boundary
	elseif indices[1] == gridsize[1]
		return f2[:,indices[1]-1,indices[2]]
	
	# If outlet is on south boundary
	elseif indices[2] == 1
		return f2[:,indices[1],indices[2]+1]

	# If outlet is on north boundary
	elseif indices[2] == gridsize[2]
		return f2[:,indices[1],indices[2]-1]
	end
end
