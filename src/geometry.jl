module geometry

export is_in_pipe, is_in_half_space, is_in_sphere

# Returns distance from cylinder center, squared
# query: DenseVector for which to compute cylindrical norm
# d₀: Direction of the cylinder
# center: Center of cylinder goes through this point
function cylindrical_norm_squared{T}(query::DenseVector{T}, d₀::DenseVector{T}, center::DenseVector{T})
    @assert length(d₀) == length(center) == length(query)
    @assert norm(d₀, 1) > 1e-8
    q₀ = query - center
    q⟂ = q₀ - (q₀⋅d₀)d₀ / (d₀⋅d₀)
    q⟂⋅q⟂
end

# True if query point is inside cylinder
# query: DenseVector for which to compute cylindrical norm
# d₀: Direction of the cylinder
# center: Center of cylinder goes through this point
# radius: of the cylinder
function is_in_pipe{T}(query::DenseVector{T}, d₀::DenseVector{T}, center::DenseVector{T}, radius::T)
    cylindrical_norm_squared(query, d₀, center) <= (radius^2)
end

# True if in half-space defined by normal and offset
# query: DenseVector for which to compute condition
# d₀: Normal to the plane defining the half-space
# origin: Origin of the plane defining the half-space
function is_in_half_space{T}(query::DenseVector{T}, d₀::DenseVector{T}, origin::DenseVector{T})
    @assert length(d₀) == length(origin) == length(query)
    @assert norm(d₀, 1) > 1e-8
    (query - origin)⋅d₀ >= 0
end


# True if in sphere defined by center and radius
# query: DenseVector for which to compute condition
function is_in_sphere{T}(query::DenseVector{T}, center::DenseVector{T}, radius::T)
    @assert length(center) == length(query)
    q₀ = center - query
    q₀⋅q₀ <= radius^2
end

end # module geometry
