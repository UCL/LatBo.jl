using Base: Cartesian

type SingleRelaxationTime{T <: FloatingPoint, DIMS} <: LatticeBoltzmann
    τ⁻¹::Float64  # Inverse of the relaxation time
    kernel::Module

    populations :: Array
    next_populations :: Array
    playground :: Array
end

function SingleRelaxationTime{T}(τ⁻¹::T, dimensions::(Int...))
    @assert τ⁻¹ > 0
    @assert length(dimensions) == 2 || length(dimensions) == 3

    if length(dimensions) == 2
        SingleRelaxationTime{T, 2}(
            τ⁻¹::T,
            D2Q9,
            zeros(T, tuple(9, dimensions...)),
            zeros(T, tuple(9, dimensions...)),
            zeros(playground.Feature, dimensions)
        )
    else
        SingleRelaxationTime{T, 3}(
            τ⁻¹::T,
            D3Q19,
            zeros(T, tuple(19, dimensions...)),
            zeros(T, tuple(19, dimensions...)),
            zeros(playground.Feature, dimensions)
        )
    end
end

for N = 2:3
    function_name = symbol("_lattice_loop_impl_$(N)D")
    @eval begin
        function $function_name(site_func::Function, sim::LatticeBoltzmann)
            playground = sim.playground
            Cartesian.@nloops $N i d->1:size(playground, d) begin
                indices = [(Cartesian.@ntuple $N i)...]
                sim.populations[:, indices...] = site_func(
                    indices, sim.populations[:, indices...],
                    Cartesian.@nref $N playground i
                )
            end
        end
    end
end

#= Iterates over the lattice sites
    * site_func: A function with signature f(indices, fᵢ, site_type). It is
        called for each lattice site.
        - indices is an array of indices for that lattice site
        - fᵢ is the array of particle probabilities for that site
        - site_type the type at that lattice site

    * sim A lattice boltzmann object with all the attendant data
=#
function lattice_loop(site_func::Function, sim::LatticeBoltzmann)
    if ndims(sim.playground) == 2
        _lattice_loop_impl_2D(site_func, sim)
    elseif ndims(sim.playground) == 3
        _lattice_loop_impl_3D(site_func, sim)
    else
        error("Undimensional LB")
    end
end
