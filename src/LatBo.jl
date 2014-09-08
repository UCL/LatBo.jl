module LatBo

export geometry, playground, LatticeBoltzmann, SingleRelaxationTime

abstract LatticeBoltzmann

include("geometry.jl")
include("playground.jl")
include("single_relaxation_time.jl")

end # module
