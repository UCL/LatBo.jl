using FactCheck: facts, context, @fact, not
#using LatBo.geometry: is_in_pipe, is_in_half_space
#using LatBo.playground: initialize, SOLID, FLUID, INLET, OUTLET, NOTHING
using LatBo: plot_frame
using Gadfly
using DataFrames

facts("creating frame elements")do

 # Create test pipe
 
 testpipe = ones(40,40)
 
 #call to function for testing, return DataFrame
 testframe=plot_frame(testpipe)
        
#obtain dimensions of frame elements
    nrows = size(testframe, 1)
    ncols = size(testframe, 2)
    
#test frame elements
@fact nrows => 41
@fact ncols => 3

end