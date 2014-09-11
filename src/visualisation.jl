module visualization
#declare use of Gadfly and DataFrames
using Gadfly
using DataFrames


###function to create dataframe from pipeline for plotting ###

function plot_frame(playpen)

	#create empty 1D arrays 
	x,y,z=Int8[],Int8[],Int8[]

	#create 1D vectors from 2D array playpen
	for i=1:size(playpen,1),j=1:size(playpen,2)
    	push!(x,i)
    	push!(y,j)
    	push!(z,playpen[i,j])
	end

	#create DataFrame
	df = DataFrame(A = x[:], B = y[:],C=z[:])

	return df

end

####function for 2d plot of pipeline####

function visualise_basic(playpen)

	#create DataFrame
	df = plot_frame(play_pen)

	#create simple 2D plot using discrete colours
	plot(df,x="A",y="B",color="C",Geom.point,Scale.discrete_color())

	end

end