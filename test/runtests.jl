using LatBo
using FactCheck: @runtest

@runtest LatBo geometry
@runtest LatBo playground
@runtest LatBo initialization
@runtest LatBo collision
@runtest LatBo integer_calc
@runtest LatBo thermodynamics
@runtest LatBo lattice_loop
@runtest LatBo zou_he_boundary
@runtest LatBo initial_probability
@runtest LatBo noslipboundary
#= @runtest LatBo one_step =#
#= @runtest LatBo plot_frame =#
#= @runtest LatBo plot_vectors =#
