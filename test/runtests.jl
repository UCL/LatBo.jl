module LatBoTests
using LatBo
using FactCheck: facts, context, @fact, not, roughly, exitstatus, exactly, @fact_throws

include("indices/cached.jl")
include("units.jl")
include("playground.jl")
include("geometry.jl")
include("indices/indices.jl")
include("lb/lattice.jl")
include("lb/thermodynamics.jl")
include("lb/collision.jl")
include("lb/streaming.jl")
include("lb/iolet.jl")
include("lb/initialization.jl")
include("lb/local_kernel.jl")
exitstatus()
end
