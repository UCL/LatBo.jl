type SingleRelaxationTime{T <: AbstractFloat} <: Collision
    # Inverse of the relaxation time
    τ⁻¹::T
end

# Defines collision kernels for SRT
collision(τ⁻¹, fᵢ, feq) = τ⁻¹ * (feq - fᵢ)
collision(k::SingleRelaxationTime, fᵢ, feq) = collision(k.τ⁻¹, fᵢ, feq)
collision!(fᵢ, k::Collision, feq, site) = (fᵢ[:, site] += collision(k, fᵢ[:, site], feq))
# Applies collision operator onto fᵢ
collision!(fᵢ, k::SingleRelaxationTime, feq, site) = collision!(fᵢ, k.τ⁻¹, feq, site)
function collision!(fᵢ, τ⁻¹, feq, site)
    const n = size(fᵢ, 1) * (site - 1)
    for i = 1:length(feq)
        fᵢ[n+i] += collision(τ⁻¹, fᵢ[n+i], feq[i])
    end
end
