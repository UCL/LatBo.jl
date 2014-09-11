using FactCheck: facts, context, @fact, not
#using LatBo.geometry: is_in_pipe, is_in_half_space
#using LatBo.playground: initialize, SOLID, FLUID, INLET, OUTLET, NOTHING
using LatBo.visualisation: plot_vectors
using Gadfly
using DataFrames

facts("creating datavector for use in visualise")do

 # Create test pipe
 
 testpipe = ones(40,40)
 
 #call to function for testing, return vectors
 x,y,z=plot_vectors(testpipe)
        
#obtain dimensions of frame elements
    xsize = size(x,1)
    ysize = size(y,1)
    zsize = size(z,1)
    
#test frame elements
@fact xsize => 40
@fact ysize => 40
@fact zsize => 40

end