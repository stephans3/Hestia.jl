#############################################################

# Segmentation

#############################################################

abstract type  AbstractSegmentation end



#############################################################

# Simple Segmentation

#############################################################


"""
    SimpleSegment
A container that identifies a list of cells with its thermal property.
"""
mutable struct SimpleSegment <: AbstractSegmentation 
    property :: AbstractMaterialProperty 
    cellindices  :: Array{T, 1} where T <: Integer

    function SimpleSegment( property :: AbstractMaterialProperty, cellindices :: Vector{<:Integer})
        if length(cellindices) == 0
            error("Array cellindices is empty! It needs the indices of heat cells.")
        end
        new(property, cellindices)
    end
end


function SimpleSegment(property :: AbstractMaterialProperty, geometry :: HeatRod)
    Nx = geometry.heatcells[1];
    cellindices     = [_ for _=1:Nx]
    return SimpleSegment(property, cellindices)
end


function SimpleSegment(property :: AbstractMaterialProperty, geometry :: HeatPlate)
    Nx = geometry.heatcells[1];
    Ny = geometry.heatcells[2];
    cellindices     = [_ for _=1:Nx*Ny]
    return SimpleSegment(property, cellindices)
end



function SimpleSegment(property :: AbstractMaterialProperty, geometry :: HeatCuboid)
    Nx = geometry.heatcells[1];
    Ny = geometry.heatcells[2];
    Nz = geometry.heatcells[3];
    cellindices     = [_ for _=1:Nx*Ny*Nz]
    return SimpleSegment(property, cellindices)
end




"""
    createSimpleSegment( property :: AbstractIsotropicProperty, cellindices :: Array{S2, 1} where S2 <: Integer)

Returns a SimpleSegment
"""
function createSimpleSegment( property :: AbstractIsotropicProperty, cellindices :: Vector{S2} where S2 <: Integer)
    if length(cellindices) == 0
        error("Array cellindices is empty! It needs the indices of heat cells.")
    end
    return SimpleSegment(property, cellindices)
end


function createSimpleSegment( property :: AbstractAnisotropicProperty, cellindices :: Vector{S2} where S2 <: Integer)
    if length(cellindices) == 0
        error("Array cellindices is empty! It needs the indices of heat cells.")
    end
    return SimpleSegment(property, cellindices)
end





# Not yet implemented

"""
    MixedSegment <: AbstractSegmentation 

A container that identifies each element of a list of cells with its own specific thermal property.
"""
mutable struct MixedSegment <: AbstractSegmentation 
    heatProperties :: Array{T1,1} where T1 <: AbstractMaterialProperty
    cellIndices    :: Array{T2,1} where T2 <: Integer
end

"""
    HyperSegment <: AbstractSegmentation 

A meta container segment to store all segments
"""
mutable struct HyperSegment <: AbstractSegmentation 
    segmentations :: Array{T,1} where T <: AbstractSegmentation
end