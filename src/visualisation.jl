module visualisation
#declare use of Gadfly, DataFrames and Grid
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

	#create DataFrame and return
	df = DataFrame(A = x[:], B = y[:],C=z[:])
	return df

end

####function for 2d plot of pipeline####

function plot2d(playpen;xlabel="X",ylabel="Y",datalabel="Type")

	#create DataFrame
	df = plot_frame(playpen)

	#create simple 2D plot using discrete colours
	plot(df,x="A",y="B",color="C",Geom.point,Scale.discrete_color(),Guide.xlabel(xlabel),Guide.ylabel(ylabel),Guide.colorkey(datalabel))

end


###function for contour plot of pipe###
function plotcont(playpen)

#create contour plot
	plot(z=playpen,Geom.contour(levels=4))

end

end #of module visualisation