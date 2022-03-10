abstract type AbstractConfiguration end
abstract type AbstractStaticConfiguration <: AbstractConfiguration end

abstract type AbstractCharacterization end


# m * exp( -||M (x - xₛ)||^ν )
"""
    RadialConfiguration <: AbstractStaticConfiguration

Stores the values for the calculation of ``m ~ exp( -||M (x - xₛ)||^{2ν} )``

`scaling` : m ∈ [0, 1]

`power` : ν ∈ [0, ∞)
    
`center` : xₛ ∈ R^{3}
    
`curvature` : M ∈ R^{3 x 3} for the planar boundaries
"""
mutable struct RadialConfiguration <: AbstractStaticConfiguration

    scaling     :: Real                         # Scaling factor:      m  ∈ [0, 1]
    power       :: Integer                      # Power of exponent :  ν  ∈ [0, ∞)
    center      :: Tuple{Real,Real,Real}        # Central point:       xₛ ∈ R^{3}
    curvature   :: Matrix{T} where T <: Real    # Curvature matrix:    M  ∈ R^{3 x 3} for the planar boundaries
end



"""
    initConfiguration()

Initializes and returns a basic `RadialConfiguration`
"""
function initConfiguration()
    scale = 1.0
    power = 1
    center = (0.0, 0.0, 0.0)
    curvature = 0.0

    return setConfiguration(scale, power, center, curvature) 
end


"""
    setConfiguration( scale :: Real , power :: Integer, curvature :: Real)

    The center is set to the origin (0,0,0).
"""
function setConfiguration( scale :: Real , power :: Integer, curvature :: Real)

    if scale < 0 || scale > 1
        println( "Scale = $(scale) is not in interval [0,1]! Scale is set to 1!" )
        scale = eltype(scale)(1)
    end
    central_point = (0.0, 0.0, 0.0);
    curvature = curvature * [1 0 0; 0 1 0; 0 0 1];

    return RadialConfiguration(scale, power, central_point, curvature)
end


"""
    setConfiguration( scale :: Real , power :: Integer; curvature = 1.0 :: Real)

    The center is set to the origin (0,0,0).
"""
function setConfiguration( scale :: Real , power :: Integer, curvature :: Matrix{T} where T <: Real)

    if scale < 0 || scale > 1
        println( "Scale = $(scale) is not in interval [0,1]! Scale is set to 1!" )
        scale = eltype(scale)(1)
    end
    central_point = (0.0, 0.0, 0.0);

    return RadialConfiguration(scale, power, central_point, curvature)
end

"""
    setConfiguration( scale :: Real , power :: Integer, central_point :: Tuple{Real,Real,Real}, curvature :: Real)

Here the variable `curvature` is multiplied with the identity matrix to gain matrix `M`.

Returns a `RadialConfiguration`
"""
function setConfiguration( scale :: Real , power :: Integer, central_point :: Tuple{Real,Real,Real}, curvature :: Real)

    if scale < 0 || scale > 1
        println( "Scale = $(scale) is not in interval [0,1]! Scale is set to 1!" )
        scale = eltype(scale)(1)
    end

    curvature = curvature * [1 0 0; 0 1 0; 0 0 1];

    return RadialConfiguration(scale, power, central_point, curvature)
end


"""
    setConfiguration( scale :: Real , power :: Integer, central_point :: Tuple{Real,Real,Real}, curvature :: Array{T,2} where T <: Real)

Returns a `RadialConfiguration`
"""
function setConfiguration( scale :: Real , power :: Integer, central_point :: Tuple{Real,Real,Real}, curvature :: Matrix{T} where T <: Real)

    if scale < 0 || scale > 1
        println( "Scale = $(scale) is not in interval [0,1]! Scale is set to 1!" )
        scale = eltype(scale)(1)
    end

    if size(curvature) != (3,3)
        error("Variable curvature has to be a scalar or a 3x3 matrix. curvature = $(curvature)")
    end

    return RadialConfiguration(scale, power, central_point, curvature)
end



function characterize(start :: T, stop :: T, step :: T, uniscale :: S) where {T <: Real, S <: Real}
    if !(0 <= uniscale <= 1)
        error("Variable uniscale (= $(uniscale)) is not inside of interval [0, 1]!")
    end

    span = range(start, step=step, stop=stop )  
    spatial_char = uniscale*ones(length(span))

    return spatial_char
end

"""
    characterize(start :: Real, stop :: Real, step :: Real, config :: RadialConfiguration; dim = 1 :: Integer )

x₁: dim = 1

x₂: dim = 2

x₃: dim = 3
"""
function characterize(start :: Real, stop :: Real, step :: Real, config :: RadialConfiguration; dim = 1 :: Integer )
    
    if dim < 1 || dim > 3
        error("Variable dim=$(dim) has to be 1, 2 or 3!")
    end
    
    span = range(start, step=step, stop=stop )  
    
    xc  = config.center[dim]
    A   = config.scaling
    M   = config.curvature
    ν   = config.power

    spatial_char = A * exp.( -abs.(M[dim,dim] * (span .- xc)).^(2*ν) )

    return spatial_char
end

function characterize(start :: Tuple{T,T}, stop :: Tuple{T,T}, step :: Tuple{T,T}, uniscale :: S) where {T <: Real, S <: Real}
    if !(0 <= uniscale <= 1)
        error("Variable uniscale (= $(uniscale)) is not inside of interval [0, 1]!")
    end
    xspan = range(start[1], step=step[1], stop=stop[1])
    yspan = range(start[2], step=step[2], stop=stop[2])
    
    spatial_char = uniscale * ones(length(xspan), length(yspan))

    return spatial_char
end


"""
    characterize(start :: Real, stop :: Real, step :: Real, config :: RadialConfiguration; dim = 1 :: Integer )

x₁ ∪ x₂ dim = {1,2}

x₁ ∪ x₃: dim = {1,3}

x₂ ∪ x₃: dim = {2,3}
"""
function characterize(start :: Tuple{T,T}, stop :: Tuple{T,T}, step :: Tuple{T,T}, config :: RadialConfiguration; dim = (1,2) :: Tuple{Integer,Integer}) where {T <: Real}
    if dim[1]==dim[2] || dim[1]*dim[2] < 1 || dim[1] > 3 || dim[2] > 3
        error("Variable dim=$(dim) has to be {1,2}, {1,3} or {2,3}!")
    end
    
    if dim[1] > dim[2]
        dim = (dim[2], dim[1])
    end

    xspan = range(start[1], step=step[1], stop=stop[1])
    yspan = range(start[2], step=step[2], stop=stop[2])
    
    xc  = [ config.center[dim[1]]; config.center[dim[2]] ]
    A   = config.scaling
    M̃   = config.curvature
    M   = [M̃[dim[1],dim[1]] M̃[dim[1],dim[2]]; M̃[dim[2],dim[1]] M̃[dim[2],dim[2]]]
    ν   = config.power

    spatial_char = zeros(T, length(xspan), length(yspan))

    for (yindex, ypos) in enumerate(yspan)
        pos = zeros(2,length(xspan))
        pos[1,:] = collect(xspan) 
        pos[2,:] .= ypos
        dist = M*(pos .- xc)
        spatial_char[:, yindex] = A * exp.( -( sqrt.( dist[1,:].^2 + dist[2,:].^2 )  ).^(2*ν) ) # Euclidean distance
        # spatial_char[:, yindex] = A * transpose(exp.( -( mapslices(x -> norm(x, 2), dist, dims=1)  ).^ν )) # Alternative with LinearAlgebra
    end

    return spatial_char
end



