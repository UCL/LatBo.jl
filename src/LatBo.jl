module LatBo

export geometry, playground, LatticeBoltzmann,
    SingleRelaxationTime, D2Q9, D3Q19, thermodynamics

abstract LatticeBoltzmann

include("geometry.jl")
include("playground.jl")
include("single_relaxation_time.jl")
include("kernel.jl")
include("thermodynamics.jl")

end # module
