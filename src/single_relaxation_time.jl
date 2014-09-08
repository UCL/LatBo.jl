type SingleRelaxationTime{T <: FloatingPoint, DIMS} <: LatticeBoltzmann
    Δt::Float64 # Time-step
    Δx::Float64 # Lattice step
    τ⁻¹::Float64  # Inverse of the relaxation time
    kernel::Symbol

    populations :: Array
    playground :: Array
end

function SingleRelaxationTime{T}(Δt::T, Δx::T, τ⁻¹::T,
    kernel::Symbol, dimensions::(Int...))

    spatial_dims, neighbors = match(
        r"D(2|3)Q(\d+)", string(kernel)
    ).captures |> int

    @assert Δx > 0
    @assert Δt > 0
    @assert τ⁻¹ > 0
    @assert spatial_dims == 2 || spatial_dims == 3
    @assert spatial_dims == length(dimensions)

    if spatial_dims == 2
        SingleRelaxationTime{T, 2}(
            Δt::T, Δx::T, τ⁻¹::T, kernel,
            zeros(T, tuple(dimensions..., neighbors)),
            zeros(playground.Feature, dimensions)
        )
    else
        SingleRelaxationTime{T, 3}(
            Δt::T, Δx::T, τ⁻¹::T, kernel,
            zeros(T, tuple(dimensions..., neighbors)),
            zeros(playground.Feature, dimensions)
        )
    end
end

