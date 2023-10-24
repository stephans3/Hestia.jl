abstract type AbstractMaterialProperty end
abstract type AbstractIsotropicProperty <: AbstractMaterialProperty end 
abstract type AbstractAnisotropicProperty <: AbstractMaterialProperty end


#############################################################

# HeatProperty

# TODO: treat dynamic and anisotropic cases

#############################################################


mutable struct StaticIsotropic <: AbstractIsotropicProperty
    λ :: Real # Thermal conductivity
    ρ :: Real # Mass density
    c :: Real # Specific heat capacity
end


"""
    getdiffusivity(prop :: StaticIsotropic)

Returns the diffusivity of a StaticIsotropic

`α = λ/(ρ ⋅ c)`
"""
function getdiffusivity(prop :: StaticIsotropic)
    return prop.λ / (prop.ρ * prop.c)
end


"""
    DynamicIsotropic

    Temperature-dependent isotropic properties

A record for static isotropic properties.
    
    ### Elements  
    `λᵢ` : thermal conductivity coefficients 
    
    `ρᵢ` : volumetric mass density coefficients
    
    `cᵢ` : specific heat capacity coefficients
"""
mutable struct DynamicIsotropic <: AbstractIsotropicProperty
    λ :: Vector{Real} # Thermal conductivity
    ρ :: Vector{Real} # Mass density
    c :: Vector{Real} # Specific heat capacity
end

"""
    specifyproperty(θ :: Real, c :: Vector{<: Real})

    - Temperature: `θ`
    - Coefficients: `c`

    Computes     
"""

# \sum_{n=1}^{N} c_{n} \theta^{n-1}
function specifyproperty(θ :: Real, c :: Vector{<: Real})
    N = length(c)
    return mapreduce(n -> c[n]*θ^(n-1), + ,1:N)
end


"""
Static anisotropic properties
λ = diag(λx, λy, λz) = [λx 0 0; 0 λy 0 ; 0 0 λz]
"""
mutable struct StaticAnisotropic <: AbstractAnisotropicProperty
    λx :: Real # Thermal conductivity in x-direction
    λy :: Real # Thermal conductivity in y-direction
    λz :: Real # Thermal conductivity in z-direction
    ρ :: Real # Mass density
    c :: Real # Specific heat capacity
end


function StaticAnisotropic(λx, λy, ρ :: Real, c :: Real)
    return StaticAnisotropic(λx, λy, 0, ρ, c)
end

# Check if length(λ_diag)==1 or length(λ_diag) > 3
function StaticAnisotropic(λ_diag :: Vector{<:Real}, ρ :: Real, c :: Real)
    return StaticAnisotropic(λ_diag..., ρ, c)
end




mutable struct DynamicAnisotropic <: AbstractAnisotropicProperty
    λx :: Vector{Real} # Thermal conductivity in x-direction
    λy :: Vector{Real} # Thermal conductivity in y-direction
    λz :: Vector{Real} # Thermal conductivity in z-direction
    ρ :: Vector{Real} # Mass density
    c :: Vector{Real} # Specific heat capacity
end



"""
createDynamicAnisotropic(conductivity_x :: Vector{T}, conductivity_y :: Vector{T}, conductivity_z :: Vector{T},  density :: Vector{T}, capacity :: Vector{T}) where T <: Real

Returns a DynamicAnisotropic    
"""
function DynamicAnisotropic(λx :: Vector{<:Real}, 
    λy :: Vector{<:Real},                                 
    ρ :: Vector{<:Real},
    c :: Vector{<:Real})
    return DynamicAnisotropic(λx, λy, [0], ρ, c)
end









#### OLD CODE BELOW





