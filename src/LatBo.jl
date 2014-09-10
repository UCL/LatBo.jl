module LatBo

export geometry, playground, LatticeBoltzmann, SingleRelaxationTime, collision

abstract LatticeBoltzmann

include("geometry.jl")
include("playground.jl")
include("single_relaxation_time.jl")
include("collision.jl")

end # module
