abstract type AbstractPhysicalProperty end
abstract type AbstractIsotropicProperty <: AbstractPhysicalProperty end 
abstract type AbstractAnisotropicProperty <: AbstractPhysicalProperty end


#############################################################

# HeatProperty

# TODO: treat dynamic and anisotropic cases

#############################################################


"""
    StaticIsoProperty

    Temperature-independent isotropic properties

A record for static isotropic properties.
    
### Elements  
`λ` : thermal conductivity 

`ρ` : volumetric mass density 

`c` : specific heat capacity
"""
mutable struct StaticIsoProperty <: AbstractIsotropicProperty
    """
    λ - denotes the thermal conductivity
    """
    λ :: Real # Thermal conductivity
    
    """
    ρ - denotes the volumetric mass density
    """
    ρ :: Real # Mass density
    
    """
    c - denotes the specific heat capacity
    """
    c :: Real # Specific heat capacity
end


"""
    createStaticIsoProperty(conductivity :: Real, density :: Real, capacity :: Real)

Returns a StaticIsoProperty
"""
function createStaticIsoProperty(conductivity :: Real, density :: Real, capacity :: Real)
    λ = conductivity # Thermal conductivity
    ρ = density      # Mass density
    c = capacity     # Specific heat capacity

    return StaticIsoProperty(λ, ρ, c)
end

"""
    getdiffusivity(prop :: StaticIsoProperty)

Returns the diffusivity of a StaticIsoProperty

`α = λ/(ρ ⋅ c)`
"""
function getdiffusivity(prop :: StaticIsoProperty)
    return prop.λ / (prop.ρ * prop.c)
end



"""
    DynamicIsoProperty

    Temperature-dependent isotropic properties

A record for static isotropic properties.
    
    ### Elements  
    `λᵢ` : thermal conductivity coefficients 
    
    `ρᵢ` : volumetric mass density coefficients
    
    `cᵢ` : specific heat capacity coefficients
"""
mutable struct DynamicIsoProperty <: AbstractIsotropicProperty
    λ :: Array{T,1} where T <: Real # Thermal conductivity
    ρ :: Array{T,1} where T <: Real # Mass density
    c :: Array{T,1} where T <: Real # Specific heat capacity
end

"""
    createDynamicIsoProperty(conductivity :: Array{T,1}, density :: Array{T,1}, capacity :: Array{T,1})  where T <: Real

Returns a DynamicIsoProperty    
"""
function createDynamicIsoProperty(conductivity :: Vector{T}, density :: Vector{T}, capacity :: Vector{T}) where T <: Real
    λ = conductivity # Thermal conductivity
    ρ = density      # Mass density
    c = capacity     # Specific heat capacity

    return DynamicIsoProperty(λ, ρ, c)
end

function specifyproperty!(prop :: Vector{T1}, temperatures :: AbstractArray{T2,1}, coefficients :: Vector{T3})  where {T1 <: Real, T2 <: Real, T3 <: Real}
    dim_temp = length(temperatures)
    dim_coeff = length(coefficients)
    tempexpo = zeros(dim_temp, dim_coeff)

    tempexpo[:,1] = ones(T2, dim_temp)

    for i = 1 : dim_coeff-1
        tempexpo[:,i+1] = temperatures.^i
    end

    LinearAlgebra.mul!(prop, tempexpo, coefficients)

end



#=
"""
Static anisotropic properties
λ = [λ₁ 0 0; 0 λ₂ 0 ; 0 0 λ₃]
"""
mutable struct StaticAnisoProperty{T1 <: Array{S1,1} where S1 <: Real, T2 <: Real} <: AbstractAnisotropicProperty
    λ :: T1 # Thermal conductivity
    ρ :: T2 # Mass density
    c :: T2 # Specific heat capacity
end

"""
Temperature-dependent anisotropic properties
Dynamic anisotropic properties
λ(θ) = [λ₁(θ) 0 0; 0 λ₂(θ) 0 ; 0 0 λ₃(θ)]
"""
mutable struct DynamicAnisoProperty{T1 <: Array{S1,2} where S1 <: Real, T2 <: Array{S2,1} where S2 <: Real} <: AbstractAnisotropicProperty
    λ :: T1 # Thermal conductivity
    ρ :: T2 # Mass density
    c :: T2 # Specific heat capacity
end
=#







################# DynamicIsoProperty: Probably not used

# The physical properties (e.g. λ, ρ, c) are approached as a power series.
function specifyproperty(temperature :: Real, coeff_vector :: Array{T,1} where T <: Real)
    prop = 0;
    for i = 1 : length(coeff_vector)
        prop += coeff_vector[i]*temperature^(i-1) # Power series approach
    end
    return prop
end


# Replace by in-place function
function specifyproperty(temperature :: T where T <: Real, property :: DynamicIsoProperty)
    λ_coeff = property.λ
    ρ_coeff = property.ρ
    c_coeff = property.c

    λ_value = specifyproperty(temperature, λ_coeff)
    ρ_value = specifyproperty(temperature, ρ_coeff)
    c_value = specifyproperty(temperature, c_coeff)

    return λ_value, ρ_value, c_value
end
