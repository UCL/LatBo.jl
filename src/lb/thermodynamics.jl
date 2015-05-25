# Density at a given lattice site
density(fᵢ::DenseVector) = sum(fᵢ)
function density(fᵢ::DenseArray)
    shape = size(fᵢ)[2:end]
    reshape([sum(fᵢ[:, i]) for i = 1:prod(shape)], shape)
end
density(sim::Simulation) = density(sim.populations)
# Momentum at given site
function momentum!(μ_out::DenseVector, fᵢ::DenseVector, cᵢ::DenseMatrix)
    @assert length(μ_out) == size(cᵢ, 1) && length(fᵢ) == size(cᵢ, 2)
    μ_out[:] = 0
    for i in 1:length(μ_out)
        for j in 1:length(fᵢ)
            if cᵢ[i, j] == 1
                μ_out[i] += fᵢ[j]
            elseif cᵢ[i, j] == -1
                μ_out[i] -= fᵢ[j]
            end
        end
    end
    μ_out
end
momentum(fᵢ::DenseVector, cᵢ::DenseMatrix) = cᵢ * fᵢ
# Velocities at a given lattice site
velocity(μ::DenseVector, ρ::Number) = μ / ρ
velocity(fᵢ::DenseVector, cᵢ::DenseMatrix, ρ::Number) = velocity(momentum(fᵢ, cᵢ), ρ)
velocity(fᵢ::DenseVector, cᵢ::DenseMatrix) = velocity(fᵢ, cᵢ, density(fᵢ))
function velocity!(ν_out::DenseVector, μ::DenseVector, ρ::Number)
    ν_out[:] = μ
    scale!(ν_out, 1 / ρ)
end

# Momentum and velocity functions that take the full grid, or a simulation object
for name in [:momentum, :velocity]
    @eval begin
        function $name(fᵢ::DenseArray, cᵢ::DenseMatrix)
            const shape = tuple(size(cᵢ, 1), size(fᵢ)[2:end]...)
            result = zeros(eltype(fᵢ), shape)
            for i = 1:prod(shape[2:end])
                result[1:shape[1], i] = $name(fᵢ[:, i], cᵢ)
            end
            result
        end
        $name(sim::Simulation) = $name(sim.populations, sim.lattice.celerities)
    end
end

#= Computes the equilibrium particle distributions $f^{eq}$:

    momentum: Macroscopic momentum μ at current lattice site
    celerities: d by n matrix of celerities ē for the lattice, with d the dimensionality of the
        lattice and n the number of particle distributions.
    weights: Weights associated with each celerity
    ρ: Density

    $f^{eq}$ = weights .* [ρ + 3ē⋅μ + \frac{9}{2ρ} (ē⋅μ)² - \frac{3}{2ρ} μ⋅μ]$
=#
function equilibrium{T, I}(
        ρ::T, momentum::DenseVector{T}, celerities::DenseMatrix{I}, weights::DenseVector{T})
    # computes momentum projected on each particle celerity first
    @assert length(momentum) == size(celerities, 1)
    @assert length(weights) == size(celerities, 2)
    μ_on_ē = celerities.'momentum
    weights .* (
        ρ
        + 3μ_on_ē
        + 9/(2ρ) * (μ_on_ē .* μ_on_ē)
        - 3/(2ρ) * dot(momentum, momentum)
    )
end
equilibrium{T, I}(lattice::Lattice{T, I}, fᵢ::DenseVector{T}) =
    equilibrium(density(fᵢ), momentum(fᵢ, lattice.celerities), lattice.celerities, lattice.weights)

equilibrium{T}(lattice::Lattice, ρ::T, momentum::DenseVector{T}) =
    equilibrium(ρ, momentum, lattice.celerities, lattice.weights)
equilibrium{T}(lattice::Symbol, ρ::T, momentum::DenseVector{T}) =
    equilibrium(getfield(LB, lattice), ρ, momentum)
# Computes equilibrium without memory-allocation
# Incoming must have the right size
function equilibrium!{T, I}(
        feq::DenseVector{T}, ρ::T, momentum::DenseVector{T}, celerities::DenseMatrix{I},
        weights::DenseVector{T})
    @assert length(feq) == size(celerities, 2)
    @assert length(momentum) == size(celerities, 1)
    @assert length(weights) == size(celerities, 2)
    const inv_tworho = 1 / (2ρ)
    # Interestingly enough, celerities[:, i] is fairly slower than celerities[j:j+d], where j is
    # incremented at each step of the loop
    const inc = length(momentum)
    const d = inc - 1
    j = 1
    for i in 1:length(feq)
        const μ_on_ē = 3dot(celerities[j:j+d], momentum)
        j += inc
        feq[i] = weights[i] * (
            ρ + μ_on_ē + (μ_on_ē * μ_on_ē - 3dot(momentum, momentum)) * inv_tworho)
    end
end
