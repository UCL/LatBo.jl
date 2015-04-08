type SingleRelaxationTime{T <: FloatingPoint} <: Collision
    # Inverse of the relaxation time
    τ⁻¹::T
end

# Defines collision kernels for SRT
collision{T}(τ⁻¹::T, fᵢ::Vector{T}, feq::Vector{T}) = τ⁻¹ * (feq - fᵢ)
function collision{T}(k::SingleRelaxationTime{T}, fᵢ::Vector{T}, feq::Vector{T})
    collision(k.τ⁻¹, feq, fᵢ)
end
# Applies collision operator onto fᵢ
collision!(k::Collision, fᵢ::Vector, feq::Vector) = (fᵢ += collision(k, fᵢ, feq))
