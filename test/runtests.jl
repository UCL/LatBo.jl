module LatBoTests
using LatBo
using FactCheck

include("thermodynamics.jl")
include("kernel.jl")
FactCheck.exitstatus()
#= @runtest LatBo geometry =#
#= @runtest LatBo playground =#
#= @runtest LatBo initialization =#
#= @runtest LatBo collision =#
#= @runtest LatBo integer_calc =#
#= @runtest LatBo lattice_loop =#
#= @runtest LatBo zou_he_boundary =#
#= @runtest LatBo initial_probability =#
#= @runtest LatBo noslipboundary =#
#= @runtest LatBo one_step =#
#= @runtest LatBo plot_frame =#
#= @runtest LatBo plot_vectors =#
end
