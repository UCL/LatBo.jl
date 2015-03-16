module playground
using Base: Cartesian
using ..Feature

# Enumerates possible features on the grid:
# - FLUID: Site on which to perform standard lattice-boltzmann
# - SOLID: Wall or boundary of some kind
# - INLET: input into the simulation domain
# - OUTLET: output out of the simulation domain
# - NOTHING: convenience value when constructing grid
# - NULL: means streaming to /dev/null
const NOTHING = convert(Feature, 1)
const FLUID   = convert(Feature, 2)
const SOLID   = convert(Feature, 3)
const INLET   = convert(Feature, 4)
const OUTLET  = convert(Feature, 5)
const NULL    = convert(Feature, 6)
const NUMBER_OF_FEATURE_TYPES = convert(Feature, 6)

for N = 2:3
    @eval begin
        function initialize(characterize::Function, grid::Array{Feature, $N})
            Cartesian.@nloops $N i grid begin
                site::Feature = characterize((Cartesian.@ntuple $N i)...)
                if site != NOTHING
                    (Cartesian.@nref $N grid i) = site
                end
            end
            grid
        end
    end
end

# Goes over the playground grid, using the function to set it up
# The function is passed the indices to the playground grid for each site in
# the grid.
# It should return an instance of Feature
# If the return is equal to NOTHING, then that site is left untouched.
# Otherwise, the site is set to that value.
# This mechanism makes it possible to call initialize multiple times with
# different functions.
initialize(characterize::Function, dimensions) = initialize(
    characterize,
    zeros(Feature, dimensions)
)

end
