#############################################################

# Boundary

#############################################################



abstract type  AbstractBoundary end
abstract type  AbstractCubicBoundary <: AbstractBoundary end

# Change boundaries as in iosetup.jl: Array -> Dict

mutable struct CubicBoundary <: AbstractBoundary
    indices :: Dict{Symbol, Vector{Int64}}       # Cell indices
    emissions :: Dict{Symbol, Vector{Emission}}  # Emission: heat transfer + heat radiation
end

function initBoundary(geometry :: AbstractGeometricalObject)
    
    position_symbols = getboundarypositions( geometry )
    
    indices = Dict{Symbol, Vector{Int64}}()
    emissions = Dict{Symbol, Vector{Emission}}()

    for pos in position_symbols
        indices[pos] = getindices(geometry, cellPosition=pos)

        N = length(indices[pos])
        zeroEmission = createEmission(0, 0, 0)
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


# The following lines will be removed in future
mutable struct HeatRodBoundary <: AbstractCubicBoundary
    indices_west :: Array{Integer,1}    # Contains 1 element
    indices_east :: Array{Integer,1}

    emissions_west :: Array{Emission,1} # Contains 1 element
    emissions_east :: Array{Emission,1}
end

mutable struct HeatPlateBoundary <: AbstractCubicBoundary
    indices_west :: Array{Integer,1}
    indices_east :: Array{Integer,1}
    indices_south :: Array{Integer,1}
    indices_north :: Array{Integer,1}

    emissions_west :: Array{Emission,1}
    emissions_east :: Array{Emission,1}
    emissions_south :: Array{Emission,1}
    emissions_north :: Array{Emission,1}
end


mutable struct HeatCuboidBoundary <: AbstractCubicBoundary
    indices_west :: Array{Integer,1}
    indices_east :: Array{Integer,1}
    indices_south :: Array{Integer,1}
    indices_north :: Array{Integer,1}
    indices_underside :: Array{Integer,1}
    indices_topside :: Array{Integer,1}

    emissions_west :: Array{Emission,1}
    emissions_east :: Array{Emission,1}
    emissions_south :: Array{Emission,1}
    emissions_north :: Array{Emission,1}
    emissions_underside :: Array{Emission,1}
    emissions_topside :: Array{Emission,1}
end

# Deprecated
function initBoundary_OLD(geometry :: AbstractGeometricalObject)
    

    if isa( geometry, HeatRod )
        iwest = [1]
        ieast = [geometry.heatcells] 

        ewest = [createEmission(0, 0, 0)]
        eeast = [createEmission(0, 0, 0)]

        indices = (iwest, ieast)
        emissions = (ewest,eeast)

        return HeatRodBoundary(indices..., emissions...)

    elseif isa( geometry, HeatPlate )
        
        iwest = getindices(geometry; cellPosition = :west)
        ieast = getindices(geometry; cellPosition = :east) 
        isouth = getindices(geometry; cellPosition = :south)
        inorth = getindices(geometry; cellPosition = :north)

        zeroEmission = createEmission(0, 0, 0)
        ewest = Array{Emission}(undef,length(iwest))
        eeast = Array{Emission}(undef,length(ieast))
        esouth = Array{Emission}(undef,length(isouth))
        enorth = Array{Emission}(undef,length(inorth))

        fill!(ewest, zeroEmission)
        fill!(eeast, zeroEmission)
        fill!(esouth, zeroEmission)
        fill!(enorth, zeroEmission)
        
        indices = (iwest, ieast, isouth, inorth)
        emissions = (ewest, eeast, esouth, enorth)

        return HeatPlateBoundary(indices..., emissions...)

    elseif isa( geometry, HeatCuboid )
          
        iwest = getindices(geometry; cellPosition = :west)
        ieast = getindices(geometry; cellPosition = :east) 
        isouth = getindices(geometry; cellPosition = :south)
        inorth = getindices(geometry; cellPosition = :north)
        iunderside = getindices(geometry; cellPosition = :underside)
        itopside = getindices(geometry; cellPosition = :topside)

        zeroEmission = createEmission(0, 0, 0)
        ewest = Array{Emission}(undef,length(iwest))
        eeast = Array{Emission}(undef,length(ieast))
        esouth = Array{Emission}(undef,length(isouth))
        enorth = Array{Emission}(undef,length(inorth))
        eunderside = Array{Emission}(undef,length(iunderside))
        etopside = Array{Emission}(undef,length(itopside))

        fill!(ewest, zeroEmission)
        fill!(eeast, zeroEmission)
        fill!(esouth, zeroEmission)
        fill!(enorth, zeroEmission)
        fill!(eunderside, zeroEmission)
        fill!(etopside, zeroEmission)
        
        indices = (iwest, ieast, isouth, inorth, iunderside, itopside)
        emissions = (ewest, eeast, esouth, enorth, eunderside, etopside)

        return HeatCuboidBoundary(indices..., emissions...)

    end

end

# Deprecated
function setEmission!_OLD(boundary :: AbstractCubicBoundary, emission :: Emission, orientation :: Symbol)
    if orientation == :west
        emission_array = Array{Emission}(undef,length(boundary.emissions_west))
        fill!(emission_array, emission)

        boundary.emissions_west = emission_array

    elseif orientation == :east
        emission_array = Array{Emission}(undef,length(boundary.emissions_east))
        fill!(emission_array, emission)

        boundary.emissions_east = emission_array


    elseif orientation == :south
        emission_array = Array{Emission}(undef,length(boundary.emissions_south))
        fill!(emission_array, emission)

        boundary.emissions_south = emission_array

    elseif orientation == :north
        emission_array = Array{Emission}(undef,length(boundary.emissions_north))
        fill!(emission_array, emission)

        boundary.emissions_north = emission_array

    elseif orientation == :topside
        emission_array = Array{Emission}(undef,length(boundary.emissions_topside))
        fill!(emission_array, emission)

        boundary.emissions_topside = emission_array

    elseif orientation == :underside
        emission_array = Array{Emission}(undef,length(boundary.emissions_underside))
        fill!(emission_array, emission)

        boundary.emissions_underside = emission_array
    end
    return nothing
end

# Deprecated
function getEmission_OLD(boundary :: AbstractCubicBoundary, cellindex :: Integer, orientation :: Symbol)

    if orientation == :west
        emission_index = findfirst(isequal(cellindex), boundary.indices_west)
        emission = boundary.emissions_west[emission_index]

    elseif orientation == :east
        emission_index = findfirst(isequal(cellindex), boundary.indices_east)
        emission = boundary.emissions_east[emission_index]

    elseif orientation == :south
        emission_index = findfirst(isequal(cellindex), boundary.indices_south)
        emission = boundary.emissions_south[emission_index]

    elseif orientation == :north
        emission_index = findfirst(isequal(cellindex), boundary.indices_north)
        emission = boundary.emissions_north[emission_index]

    elseif orientation == :topside
        emission_index = findfirst(isequal(cellindex), boundary.indices_topside)
        emission = boundary.emissions_topside[emission_index]

    elseif orientation == :underside
        emission_index = findfirst(isequal(cellindex), boundary.indices_underside)
        emission = boundary.emissions_underside[emission_index] 

    end

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
