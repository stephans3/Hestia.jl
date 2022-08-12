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
    heatProperty :: AbstractPhysicalProperty 
    cellIndices  :: Array{T, 1} where T <: Integer
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

"""
    MixedSegment <: AbstractSegmentation 

A container that identifies each element of a list of cells with its own specific thermal property.
"""
mutable struct MixedSegment <: AbstractSegmentation 
    heatProperties :: Array{T1,1} where T1 <: AbstractPhysicalProperty
    cellIndices    :: Array{T2,1} where T2 <: Integer
end

"""
    HyperSegment <: AbstractSegmentation 

A meta container segment to store all segments
"""
mutable struct HyperSegment <: AbstractSegmentation 
    segmentations :: Array{T,1} where T <: AbstractSegmentation
end