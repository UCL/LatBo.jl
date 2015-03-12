type Lattice{T <: Real}
    celerities::Matrix{T}
    weights::Vector{T}
end
const D2Q9 = Lattice{Float64}(
  [0 1 0 -1  0 1  1 -1 -1;
   0 0 1  0 -1 1 -1  1 -1],
  vcat(4./9., [1./9. for i=1:4], [1./36. for i=1:4])
)

const D3Q19 = Lattice{Float64}(
  [0 1 0 0 -1  0  0 1 -1  1 -1 1 -1  1 -1 0  0  0  0
   0 0 1 0  0 -1  0 1 -1 -1  1 0  0  0  0 1 -1  1 -1
   0 0 0 1  0  0 -1 0  0  0  0 1 -1 -1  1 1 -1 -1  1],
  vcat(1./3., [1./18. for i=1:6], [1./36. for i=7:18])
)
