module geometry

export is_in_pipe, is_in_half_space, is_in_sphere

# Returns distance from cylinder center, squared
# query: Vector for which to compute cylindrical norm
# d₀: Direction of the cylinder
# center: Center of cylinder goes through this point
function cylindrical_norm_squared{T}(query::Array{T, 1}, d₀::Array{T, 1},
    center::Array{T, 1})
    @assert length(d₀) == length(center) == length(query)
    @assert norm(d₀, 1) > 1e-8
    q₀ = query - center
    q⟂ = q₀ - (q₀⋅d₀)d₀ / (d₀⋅d₀)
    q⟂⋅q⟂
end

# True if query point is inside cylinder
# query: Vector for which to compute cylindrical norm
# d₀: Direction of the cylinder
# center: Center of cylinder goes through this point
# radius: of the cylinder
function is_in_pipe{T}(query::Array{T, 1}, d₀::Array{T, 1},
    center::Array{T, 1}, radius::T)
    cylindrical_norm_squared(query, d₀, center) <= (radius^2)
end

# True if in half-space defined by normal and offset
# query: Vector for which to compute condition
# d₀: Normal to the plane defining the half-space
# origin: Origin of the plane defining the half-space
function is_in_half_space{T}(query::Array{T, 1}, d₀::Array{T, 1},
    origin::Array{T, 1})
    @assert length(d₀) == length(origin) == length(query)
    @assert norm(d₀, 1) > 1e-8
    (query - origin)⋅d₀ >= 0
end


# True if in sphere defined by center and radius
# query: Vector for which to compute condition
function is_in_sphere{T}(query::Array{T, 1}, center::Array{T, 1}, radius::T)
    @assert length(center) == length(query)
    q₀ = center - query
    q₀⋅q₀ <= radius^2
end

end # module geometry
