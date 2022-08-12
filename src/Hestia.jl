# Author: Stephan Scholz
# Year 2022

module Hestia

using Base: Tuple, Integer


# Core utilities
export StaticIsoProperty, DynamicIsoProperty
export createStaticIsoProperty, createDynamicIsoProperty, getdiffusivity
export StaticAnisoProperty, DynamicAnisoProperty
export createStaticAnisoProperty, createDynamicAnisoProperty
include("core/properties.jl")

export AbstractSegmentation, SimpleSegment, MixedSegment, HyperSegment
export createSimpleSegment
include("core/segmentation.jl")

export HeatRod, HeatPlate, HeatCuboid
export getHeatCells, getindices
#export getnumofboundarycells, getboundarypositions, getpositionofindex, getindexofposition # Only for testing
include("core/geometry.jl")


export Emission
export createEmission, emit
include("core/emission.jl")


#export HeatRodBoundary, HeatPlateBoundary, HeatCuboidBoundary
export HeatCuboidBoundary
export initBoundary, setEmission!, getEmission
include("core/boundary.jl")


#### Interface #####
export RadialConfiguration
export initConfiguration, setConfiguration, characterize
include("interface/characterization.jl")

export IOSetup
export initIOSetup, setIOSetup!, getActuation, getSensing, measureWAM
export checkIOSetup2D
include("interface/iosetup.jl")

# 
#export induce!
export HeatRodActuation, HeatPlateActuation, HeatCuboidActuation
export initActuation, getActuation, setActuation!, checkActuation2D
include("interface/actuators.jl")

# Core
export CubicHeatProblem
export diffusion!
include("core/heatproblem.jl")




end # module
