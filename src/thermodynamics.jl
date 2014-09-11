module thermodynamics

# Density at a given lattice site
density{T}(fᵢ::Array{T, 1}) = sum(fᵢ)
# Velocities at a given lattice site
velocity{T}(fᵢ::Array{T, 1}, cᵢ::Array{T, 2}, ρ::T) = (
    vec(sum(cᵢ .* transpose(fᵢ), 2) / ρ)
)
velocity{T}(fᵢ::Array{T, 1}, cᵢ::Array{T, 2}) = velocity(fᵢ, cᵢ, density(fᵢ))
# deviatoric tensor at given point
function deviatoric{T}(fᵢ::Array{T, 1}, fᵢ⁼::Array{T, 1}, cᵢ::Array{T, 2},
    τ⁻¹::T)
   σ = zeros(Float64, size(cᵢ, 1), size(cᵢ, 1))
   for j = 1:size(cᵢ, 2)
       σ += (fᵢ[j] - fᵢ⁼[j]) .* (cᵢ[:, j] .* transpose(cᵢ[:, j]))
   end
   (1. - 0.5τ⁻¹)σ
end

end
