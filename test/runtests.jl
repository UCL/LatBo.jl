module LatBoTests
using FactCheck: facts, context, @fact, not, roughly, exitstatus, exactly

include("geometry.jl")
include("indexing.jl")
include("lb/lattice.jl")
include("lb/thermodynamics.jl")
include("lb/collision.jl")
include("lb/streaming.jl")
include("lb/iolet.jl")
include("lb/local_kernel.jl")
exitstatus()
end
