using FactCheck: facts, context, @fact, not
#using LatBo.geometry: is_in_pipe, is_in_half_space
#using LatBo.playground: initialize, SOLID, FLUID, INLET, OUTLET, NOTHING
using Gadfly
using DataFrames

facts("creating frame elements")do

 # Create test pipe
 
 testpipe = ones(40,40)
 
 #call to function for testing, return DataFrame
 testframe=visualisation_basic(testpipe)
        
#Test dimensions of frame elements
    