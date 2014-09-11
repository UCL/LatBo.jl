module LatBo

export geometry, playground, LatticeBoltzmann,
    SingleRelaxationTime, D2Q9, D3Q19, thermodynamics, collision,
	noslip_boundary

abstract LatticeBoltzmann

include("geometry.jl")
include("playground.jl")
include("single_relaxation_time.jl")
include("collision.jl")
include("integer_calc.jl")
include("kernel.jl")
include("thermodynamics.jl")
include("noslip_boundary.jl")

end # module
