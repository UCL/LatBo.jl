type SingleRelaxationTime{T <: FloatingPoint} <: Collision
    # Inverse of the relaxation time
    τ⁻¹::T
end

# Defines collision kernels for SRT
collision{T}(τ⁻¹::T, fᵢ::Vector{T}, feq::Vector{T}) = τ⁻¹ * (feq - fᵢ)
collision(k::SingleRelaxationTime, fᵢ::Vector, feq::Vector) = collision(k.τ⁻¹, fᵢ, feq)
# Applies collision operator onto fᵢ
collision!(k::Collision, fᵢ::SubArray, feq::Vector) = (fᵢ += collision(k, fᵢ, feq))
