type SandBox{T <: AbstractFloat, I <: Integer} <: Simulation{T, I}
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
function SandBox{T, N, I <: Integer}(lattice::LB.Lattice{T, I},
                                     dimensions::NTuple{N, I};
                                     kwargs...)
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
    SandBox(
        lattice::LB.Lattice{T, I},
        Indices.Cached(Indices.Cartesian(I[dimensions...]), lattice),
        getarg(:kernels, Dict{Playground.Feature, LB.LocalKernel}()),
        initializers,
        getarg(:population, zeros(T, tuple(n, dimensions...))),
        getarg(:next_population, zeros(T, tuple(n, dimensions...))),
        getarg(:playground, Playground.FLUID * ones(Playground.Feature, dimensions...)),
        getarg(:time, 0)
    )
end
function SandBox(lattice::Symbol, args::NTuple; kwargs...)
    lattice = getfield(LB, lattice)::LB.Lattice
    SandBox(lattice, args; kwargs...)
end
