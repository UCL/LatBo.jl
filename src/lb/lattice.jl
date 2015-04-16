type Lattice{T <: Real, I <: Int}
    celerities::Matrix{I}
    weights::Vector{T}
    inversion::Vector{I}

    # Creates inversion list automatically
    function Lattice(celerities::Matrix{T}, weights::Vector{T})
        @assert(size(celerities, 2) == length(weights))
        inversion = Int64[-1 for i in 1:length(weights)]
        for i in 1:length(weights)
            if inversion[i] != -1
                continue
            end
            for j in i:length(weights)
                if all(abs(celerities[:, i] + celerities[:, j]) .< 1e-8)
                   inversion[i], inversion[j] = j, i
                   break
                end
            end
        end
        @assert all(inversion .!= -1)

        new(celerities, weights, inversion)
    end
end


const D2Q9 = Lattice{Float64, Int64}(
  Float64[0 1 0 -1  0 1  1 -1 -1;
          0 0 1  0 -1 1 -1  1 -1],
  vcat(4./9., [1./9. for i=1:4], [1./36. for i=1:4])
)

const D3Q19 = Lattice{Float64, Int64}(
  Float64[0 1 0 0 -1  0  0 1 -1  1 -1 1 -1  1 -1 0  0  0  0
          0 0 1 0  0 -1  0 1 -1 -1  1 0  0  0  0 1 -1  1 -1
          0 0 0 1  0  0 -1 0  0  0  0 1 -1 -1  1 1 -1 -1  1],
  vcat(1./3., [1./18. for i=1:6], [1./36. for i=7:18])
)

const speed_of_sound_squared = 1.0 / 3.0
const speed_of_sound = sqrt(speed_of_sound_squared)

# Expand usage of neighbor_index
function neighbor_index(
        indexing::Indexing, site::Union(GridCoords, Integer), lattice::Lattice, direction::Integer)
    neighbor_index(indexing, site, lattice.celerities[:, direction])
end
function neighbor_index(sim::Simulation, site::Union(Integer, GridCoords), direction::Integer)
    neighbor_index(sim.indexing, site, sim.lattice, direction)
end

length(lattice::Lattice) = length(lattice.weights)
