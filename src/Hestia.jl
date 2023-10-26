# Author: Stephan Scholz
# Year 2022

module Hestia

using Base: Tuple, Integer


# Core utilities
export StaticIsotropic, DynamicIsotropic
export StaticAnisotropic, DynamicAnisotropic
export getdiffusivity, specifyproperty
include("core/properties.jl")


export HeatRod, HeatPlate, HeatCuboid
export getHeatCells, getindices
#export getnumofboundarycells, getboundarypositions, getpositionofindex, getindexofposition # Only for testing
include("core/geometry.jl")

# export AbstractSegmentation, SimpleSegment, MixedSegment, HyperSegment
# export createSimpleSegment
# include("core/segmentation.jl")


export Emission
export emit, emit!
include("core/emission.jl")


export HeatCuboidBoundary
export Boundary, setEmission!, getEmission
include("core/boundary.jl")


#### Interface #####
export RadialCharacteristics
export characterize
include("interface/characterization.jl")

export IOSetup
export IOSetup, setIOSetup!, getActuation, getSensing, measure, getCharacteristics
export checkIOSetup2D
include("interface/iosetup.jl")

# 
#export induce!
# export HeatRodActuation, HeatPlateActuation, HeatCuboidActuation
# export initActuation, getActuation, setActuation!, checkActuation2D
# include("interface/actuators.jl")

# Core
export CubicHeatProblem
export diffusion!
include("core/heatproblem.jl")

#export diffusion_fast!, diffusion_fast_xy!
#include("experimental/heatproblem_fast.jl")

end # module
