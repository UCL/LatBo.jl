module thermodynamics

# Density at a given lattice site
density(fᵢ) = sum(fᵢ)
# Velocities at a given lattice site
velocity(fᵢ, cᵢ, ρ) = vec(sum(cᵢ .* transpose(fᵢ), 2))# / ρ)
velocity(fᵢ, cᵢ) = velocity(fᵢ, cᵢ, density(fᵢ))
# deviatoric tensor at given point
function deviatoric(fᵢ, fᵢ⁼, cᵢ, τ⁻¹)
   σ = zeros(Float64, size(cᵢ, 1), size(cᵢ, 1))
   for j = 1:size(cᵢ, 2)
       σ += (fᵢ[j] - fᵢ⁼[j]) .* (cᵢ[:, j] .* transpose(cᵢ[:, j]))
   end
   (1. - 0.5τ⁻¹)σ
end

end
