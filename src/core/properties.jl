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
    λ :: Real # Thermal conductivity
    ρ :: Real # Mass density
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
    λ :: Vector{Real} # Thermal conductivity
    ρ :: Vector{Real} # Mass density
    c :: Vector{Real} # Specific heat capacity
end

"""
createDynamicIsoProperty(conductivity :: Vector{T1}, density :: Vector{T2}, capacity :: Vector{T3}) where where {T1 <: Real, T2 <: Real, T3 <: Real}

Returns a DynamicIsoProperty    
"""
function createDynamicIsoProperty(conductivity :: Vector{T1}, density :: Vector{T2}, capacity :: Vector{T3}) where {T1 <: Real, T2 <: Real, T3 <: Real}
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




"""
Static anisotropic properties
λ = diag(λx, λy, λz) = [λx 0 0; 0 λy 0 ; 0 0 λz]
"""
mutable struct StaticAnisoProperty <: AbstractAnisotropicProperty
    λx :: Real # Thermal conductivity in x-direction
    λy :: Real # Thermal conductivity in y-direction
    λz :: Real # Thermal conductivity in z-direction
    ρ :: Real # Mass density
    c :: Real # Specific heat capacity
end


function createStaticAnisoProperty(conductivity_x :: Real, conductivity_y :: Real, density :: Real, capacity :: Real)
    λx = conductivity_x # Thermal conductivity
    λy = conductivity_y # Thermal conductivity
    λz = 0; # Thermal conductivity

    ρ = density      # Mass density
    c = capacity     # Specific heat capacity

    return StaticAnisoProperty(λx, λy, λz, ρ, c)
end


function createStaticAnisoProperty(conductivity_x :: Real, conductivity_y :: Real, conductivity_z :: Real, density :: Real, capacity :: Real)
    λx = conductivity_x # Thermal conductivity
    λy = conductivity_y # Thermal conductivity
    λz = conductivity_z # Thermal conductivity

    ρ = density      # Mass density
    c = capacity     # Specific heat capacity

    return StaticAnisoProperty(λx, λy, λz, ρ, c)
end



mutable struct DynamicAnisoProperty <: AbstractAnisotropicProperty
    λx :: Vector{Real} # Thermal conductivity in x-direction
    λy :: Vector{Real} # Thermal conductivity in y-direction
    λz :: Vector{Real} # Thermal conductivity in z-direction
    ρ :: Vector{Real} # Mass density
    c :: Vector{Real} # Specific heat capacity
end



"""
createDynamicAnisoProperty(conductivity_x :: Vector{T}, conductivity_y :: Vector{T}, conductivity_z :: Vector{T},  density :: Vector{T}, capacity :: Vector{T}) where T <: Real

Returns a DynamicAnisoProperty    
"""
function createDynamicAnisoProperty(conductivity_x :: Vector{T1}, 
                                    conductivity_y :: Vector{T2}, 
                                    density :: Vector{T3}, 
                                    capacity :: Vector{T4}) where {T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real}
    λx = conductivity_x # Thermal conductivity in x-direction
    λy = conductivity_y # ... in y-direction
    λz = Real[]         # ... in z-direction
    ρ = density      # Mass density
    c = capacity     # Specific heat capacity

    return DynamicAnisoProperty(λx, λy, λz, ρ, c)
end

"""
createDynamicAnisoProperty(conductivity_x :: Vector{T}, conductivity_y :: Vector{T}, conductivity_z :: Vector{T},  density :: Vector{T}, capacity :: Vector{T}) where T <: Real

Returns a DynamicAnisoProperty    
"""
function createDynamicAnisoProperty(conductivity_x :: Vector{T1}, 
                                    conductivity_y :: Vector{T2}, 
                                    conductivity_z :: Vector{T3},  
                                    density :: Vector{T4}, 
                                    capacity :: Vector{T5}) where {T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real}
    λx = conductivity_x # Thermal conductivity in x-direction
    λy = conductivity_y # ... in y-direction
    λz = conductivity_z # ... in z-direction
    ρ = density      # Mass density
    c = capacity     # Specific heat capacity

    return DynamicAnisoProperty(λx, λy, λz, ρ, c)
end



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
