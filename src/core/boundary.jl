#############################################################

# Boundary

#############################################################



abstract type  AbstractBoundary end
abstract type  AbstractCubicBoundary <: AbstractBoundary end

mutable struct CubicBoundary <: AbstractBoundary
    indices :: Dict{Symbol, Vector{Int64}}       # Cell indices
    emissions :: Dict{Symbol, Vector{Emission}}  # Emission: heat transfer + heat radiation
end


"""
    Boundary(geometry :: AbstractGeometricalObject)

Initialize the boundary sides for a geometry. 
"""
function Boundary(geometry :: AbstractGeometricalObject)
    
    position_symbols = getboundarypositions( geometry )
    
    indices = Dict{Symbol, Vector{Int64}}()
    emissions = Dict{Symbol, Vector{Emission}}()

    for pos in position_symbols
        indices[pos] = getindices(geometry, cellPosition=pos)

        N = length(indices[pos])
        zeroEmission = Emission(0, 0, 0)
        emissions[pos] = Vector{Emission}(undef,N)
        fill!(emissions[pos], zeroEmission)
    end

    return CubicBoundary(indices,emissions)
end

function setEmission!(boundary :: CubicBoundary, emission :: Emission, orientation :: Symbol)

    N = length(boundary.emissions[orientation])
    emission_array = Array{Emission}(undef,N)
    fill!(emission_array, emission)

    boundary.emissions[orientation] = emission_array;

    return nothing
end

function getEmission(boundary :: CubicBoundary, cellindex :: Integer, orientation :: Symbol)

    emission_index = findfirst(isequal(cellindex), boundary.indices[orientation])
    emission = (boundary.emissions[orientation])[emission_index]

    return emission
end




# RENAME
function createBoundaryMatrix(indices, totallength)
    if maximum(indices) > totallength
        error("Total length is smaller than maximum element of array indices!")
    end

    bMatrix = spzeros(Int64, totallength, length(indices))
    
    sort!(indices)

    for (idx, val) in enumerate(indices)
        bMatrix[val, idx] = 1
    end

    return bMatrix
end
