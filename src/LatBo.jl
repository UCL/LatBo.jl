module LatBo

export SandBox, run!

# Type defining the feature of the simulation playground
abstract Simulation{T <: FloatingPoint, I <: Int}

include("indexing.jl")
include("playground.jl")
include("geometry.jl")
include("lb/lb.jl")

type SandBox{T <: FloatingPoint, I <: Int} <: Simulation{T, I}
    # Lattice on which the kernel acts
    lattice::LB.Lattice{T, I}
    # indexing kernel
    indexing::Indices.Indexing
    # Local kernels for each lattice type
    kernels :: Dict{Playground.Feature, LB.LocalKernel}
    # Initializes the populations
    initializers :: Dict{Playground.Feature, LB.Initializer}
    # Current population
    populations :: Array{T}
    # Next population
    next_populations :: Array{T}
    # Describe where flow takes place
    playground :: Array{Playground.Feature}
    # Current time step
    time::I
end

# Simple constructor for simulation structure
function SandBox{T, I}(lattice::LB.Lattice{T, I}, dimensions::(Integer...); kwargs...)
    function getarg(k::Symbol, default)
        for (key, value) in kwargs
            if key == k
                return value
            end
        end
        return default
    end

    initializers = getarg(:initializers, Dict{Playground.Feature, LB.Initializer}())
    keys = [k for (k, v) in kwargs]
    if :ρ₀ ∈ keys && :μ₀ ∈ keys
        initializers[Playground.FLUID] =
            LB.Homogeneous{T}(getarg(:ρ₀, 0), getarg(:μ₀, zeros(T, length(dimensions))))
    end

    const n = length(lattice.weights)
    SandBox{T, I}(
        lattice::LB.Lattice{T, I},
        Indices.Cartesian(I[dimensions...]),
        getarg(:kernels, Dict{Playground.Feature, LB.LocalKernel}()),
        initializers,
        getarg(:population, zeros(T, tuple(n, dimensions...))),
        getarg(:next_population, zeros(T, tuple(n, dimensions...))),
        getarg(:playground, Playground.FLUID * ones(Playground.Feature, dimensions...)),
        getarg(:time, 0)
    )
end
function SandBox(lattice::Symbol, args...; kwargs...)
    lattice = getfield(LB, lattice)::LB.Lattice
    SandBox(lattice, args...; kwargs...)
end

# Runs lattice boltzmann for single step
function run!(observer::Function, sim::Simulation; doinit::Bool=true, nsteps::Integer=1)
    # First initializes lattice
    if doinit
        LB.initialize(sim)
    end
    # Then run for N steps
    for step in 1:nsteps
        # run local kernel for each site
        LB.local_kernel(sim)
        # swap populations
        sim.populations, sim.next_populations = sim.next_populations, sim.populations
        # run observer at each step
        observer()
    end
end

# Run simulation for N numbers of steps without observing
run!(sim::Simulation; kwargs...) = run!(()->nothing, sim; kwargs...)

end # module
