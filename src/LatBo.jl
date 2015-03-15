module LatBo

export geometry, playground, LatticeBoltzmann, SingleRelaxationTime, D2Q9,
    D3Q19, thermodynamics, collision, lattice_loop, integer_calc,
    noslip_boundary, run_lb, kernel

# Base type for all kernally stuff
abstract Kernel
# Local kernel for each lattice type
abstract LocalKernel <: Kernel
# Base type for all indexing kernels
abstract Indexing <: Kernel
# Base type for all collision kernels
abstract Collision <: Kernel
# Base type for all collision kernels
abstract Streaming <: Kernel
# Type defining the feature of the simulation playground
typealias Feature Uint8

include("lattice.jl")

type Simulation{T <: Real, I <: Int}
    # Lattice on which the kernel acts
    lattice::Lattice{T, I}
    # indexing kernel
    indexing::Indexing
    # Local kernels for each lattice type
    kernels :: Dict{Feature, LocalKernel}
    # Current population
    populations :: Array{T}
    # Next population
    next_populations :: Array{T}
    # Describe where flow takes place
    playground :: Array{Feature}
    # Current time step
    time::T
end


include("geometry.jl")
include("playground.jl")
include("thermodynamics.jl")
include("kernel.jl")

# Simple constructor for simulation structure
function Simulation{T, I}(lattice::Lattice{T, I}, dimensions::(Int...);kwargs...)
    function getarg(k::Symbol, default)
        for (key, value) in kwargs
            if key == k
                return value
            end
        end
        return default
    end

    simground = getarg(:playground, zeros(Feature, dimensions...))
    if :playground ∉ kwargs
        simground[:] = playground.FLUID
    end
    const n = length(lattice.weights)
    Simulation{T, I}(
      lattice,
      lb.Cartesian(T[u for u in dimensions]),
      getarg(:kernels, Dict{Feature, LocalKernel}()),
      getarg(:population, zeros(T, tuple(n, dimensions...))),
      getarg(:next_population, zeros(T, tuple(n, dimensions...))),
      simground,
      getarg(:time, 0)
    )
end
#= include("single_relaxation_time.jl") =#
#= include("collision.jl") =#
#= include("integer_calc.jl") =#
#= include("zou_he_boundary.jl") =#
#= include("initial_probability.jl") =#
#= include("noslip_boundary.jl") =#


# Creates single-relaxation time simulation
#= function Simulation{T}( =#
#=     τ⁻¹::T, dimensions::(Int...); indexing::Type{kernel.Collision} = kernel.Cartesian) =#
#=     @assert τ⁻¹ > 0 =#
#=     @assert length(dimensions) == 2 || length(dimensions) == 3 =#
#=     @assert indexing <: kernel.IndexingKernel =#
#=  =#
#=     if length(dimensions) == 2 =#
#=         lattice = D2Q9 =#
#=     else =#
#=         lattice = D3Q19 =#
#=     end =#
#=     # Sets up the streamers =#
#=     streaming = [NullStreaming() for i in 1:playground.NUMBER_OF_FEATURE_TYPES] =#
#=     streaming[playground.FLUID] = kernel.BulkStreaming() =#
#=  =#
#=     Simulation{T}( =#
#=         lattice, =#
#=         SingleRelaxationTimeCollision(τ⁻¹), =#
#=         indexing, =#
#=         zeros(T, tuple(len(lattice.weights), dimensions...)), =#
#=         zeros(T, tuple(len(lattice.weights), dimensions...)), =#
#=         zeros(playground.Feature, dimensions) =#
#=     ) =#
#= end =#

#= # Runs lattice boltzmann for n steps =#
#= function run_lb(observer::Function, sim::LatticeBoltzmann, nsteps::Int) =#
#=     for n = 1:nsteps =#
#=         run_lb(observer, sim) =#
#=     end =#
#= end =#
#=  =#
#= # Runs lattice boltzmann for single step =#
#= function run_lb(observer::Function, sim::LatticeBoltzmann) =#
#=     # Aliases for easier acces to quantities =#
#=     this_pop = sim.populations =#
#=     next_pop = sim.next_populations =#
#=     gridsize = [size(sim.playground)...] =#
#=     kernel = sim.kernel =#
#=     celerities = sim.kernel.celerities =#
#=  =#
#=     # Loop over each lattice site =#
#=     lattice_loop(sim) do indices, fᵢ, feature =#
#=         # Apply collision step =#
#=         this_pop[:, indices...] += collision(fᵢ, kernel, sim.τ⁻¹) =#
#=         # Apply streaming step =#
#=         for v = 1:size(celerities, 2) =#
#=             streamed = integer_calc(gridsize, indices, celerities[:, v]) =#
#=             next_pop[v, streamed...] = this_pop[v, indices...] =#
#=         end =#
#=     end =#
#=  =#
#=     # Applying boundary condition =#
#=     lattice_loop(sim) do indices, fᵢ, feature =#
#=         if feature == playground.INLET =#
#=             next_pop[:, indices...] = zou_he_boundary( =#
#=                 indices..., gridsize..., fᵢ, sim.inlet_velocity) =#
#=         elseif feature == playground.SOLID =#
#=             noslip_boundary(sim.playground,indices,next_pop) =#
#=         end =#
#=     end =#
#=  =#
#=     # reassign populations and next_populations with swap =#
#=     sim.populations  = next_pop =#
#=     sim.next_populations = this_pop =#
#=  =#
#=     # run observer at each step =#
#=     observer() =#
#= end =#
#=  =#
#= # Run simulation for N numbers of steps =#
#= run_lb(sim::LatticeBoltzmann, nsteps::Int) = ( =#
#=     run_lb(()->nothing, sim, nsteps) =#
#= ) =#
#= run_lb(sim::LatticeBoltzmann) = run_lb(()->nothing, sim) =#

end # module
