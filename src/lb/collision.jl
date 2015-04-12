type SingleRelaxationTime{T <: FloatingPoint} <: Collision
    # Inverse of the relaxation time
    τ⁻¹::T
end

# Defines collision kernels for SRT
collision(τ⁻¹::Number, fᵢ::DenseVector, feq::DenseVector) = τ⁻¹ * (feq - fᵢ)
collision(k::SingleRelaxationTime, fᵢ::DenseVector, feq::DenseVector) = collision(k.τ⁻¹, fᵢ, feq)
# Applies collision operator onto fᵢ
collision!(k::Collision, fᵢ::SubArray, feq::DenseVector) = (fᵢ += collision(k, fᵢ, feq))
