# 2023 / 10 / 25
## characterization.jl
- Rename Configuration -> Characteristics
- Change constructuro initConfiguration -> RadialCharacteristics

# iosetup.jl
- Add function getCharacteristics



# 2023 / 10 / 24
## emission.jl
- Remove old constructor
- Add docs for new constructor 
- Add emit!

## boundary.jl
- change initBoundary to Boundary

## iosetup.jl
- change measureWAM to measure
- add some docs
- change initIOSetup to IOSetup

# 2023 / 10 / 23
## properties.jl
- StaticIsoProperty -> StaticIsotropic
- DynamicIsoProperty -> DynamicIsotropic
- StaticAnisoProperty -> StaticAnisotropic
- DynamicAnisoProperty -> DynamicAnisotropic
- Remove old parts
- specifyproperty now with mapreduce



## emissions.jl
- Update constructor for type Emission

## boundary.jl
- Remove deprecated parts
- change Emission constructor

## segmentation.jl
- Update constructor SimpleSegmentation
- Add constructor for HeatRod, HeatPlate, HeatCuboid

## geometry.jl
- Remove segmentation from HeatRod, HeatPlate, HeatCuboid

