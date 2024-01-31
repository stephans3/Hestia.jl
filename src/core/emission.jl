abstract type AbstractEmission end


#############################################################

# Constant heat radiation + conduction

#############################################################


"""
    Emission  <: AbstractEmission

Type Emission contains the coefficients for (convective) heat transfer and heat radiation 
    
```math
\\Phi = -h ~ (\\theta - \\theta_{amb}) - \\sigma ~ [\\varepsilon_{1} ~ \\theta^4 - \\varepsilon_{2} ~ \\theta_{amb}^4)]
```    

Constructor `Emission(h, θamb, ε₁, ε₂, θext)` expects emissivity `ε₁` and `ε₂` which must be in interval [0,1]. The Stefan-Boltzmann constant is stored internally: `σ = 5.6703744191844294e-8`.

### Elements

`h` : heat transfer coefficient

`θamb` : ambient temperature

`ε₁` : Emissivity of main object

`ε₂` : Emissivity of external second body

`θext` : temperature of an external second body
"""
mutable struct Emission <: AbstractEmission 
   h :: Real    # Heat conduction coefficient
   θamb :: Real # Ambient temperature  
   ε₁ :: Real   # Emissivity of main object
   ε₂ :: Real   # Emissivity of external second body
   θext :: Real # Temperature of an external second body

   function Emission(h, θamb, ε₁, ε₂, θext)
    if  ((0 <= ε₁ <= 1) == false)
        error("Emissivity ε₁ has to be in interval [0, 1]!")
    elseif ((0 <= ε₂ <= 1) == false) 
        error("Emissivity ε₂ has to be in interval [0, 1]!")
    elseif h < 0
        error("Heat transfer coefficient h has to be greater than or equal to zero!")
    elseif θamb < 0
        error("Ambient temperature θamb has to be greater than or equal to zero!")
    elseif θext < 0
        error("External temperature θext has to be greater than or equal to zero!")
    end

    new(h, θamb, ε₁, ε₂, θext)
   end
end

"""
    Constructor for type Emission.

Sets both heat radiation terms to zero: `ε₁=0`, `ε₂=0` and `θext=0`
```math
\\Phi = -h ~ (\\theta - \\theta_{amb}) 
```    
"""
function Emission(h, θamb)
    Emission(h, θamb, 0, 0, 0)
end 


"""
    Constructor for type Emission.

Sets external heat radiation term to zero: `ε₂=0` and `θext=0`
```math
\\Phi = -h ~ (\\theta - \\theta_{amb}) - \\sigma ~ \\varepsilon_{1} ~ \\theta^4
```    
"""
function Emission(h, θamb, ε₁)
    Emission(h, θamb, ε₁, 0, 0)
end 


"""
    emit(temperature :: Real, emission :: Emission)

Calculates the right-hand side of the boundary conditions for a given `Emission`.

```math
\\Phi = -h ~ (\\theta - \\theta_{amb}) - \\sigma ~ [\\varepsilon_{1} ~ \\theta^4 - \\varepsilon_{2} ~ \\theta_{amb}^4)]
```  
"""
function emit(temperature :: Real, emission :: Emission)
    sb = 5.6703744191844294e-8 # Stefan-Boltzmann constant

    θ = temperature
    θamb = emission.θamb
    θext = emission.θext

    h = emission.h
    ε₁ = emission.ε₁
    ε₂ = emission.ε₂

    flux = -h*(θ - θamb) - sb*(ε₁*θ^4 - ε₂*θext^4)

    return flux
end


"""
    emit!(flux :: Vector{<:Real}, temperature :: Vector{<:Real}, emission :: Emission)

Calculates the right-hand side of the natural Robin boundary along a boundary for a given `Emission`.

```math
\\Phi = -h ~ (\\theta - \\theta_{amb}) - \\sigma ~ [\\varepsilon_{1} ~ \\theta^4 - \\varepsilon_{2} ~ \\theta_{amb}^4)]
```  

Note: in-place operation - results are saved in array `flux`.
"""
function emit!(flux :: Vector{<:Real}, temperature :: Vector{<:Real}, emission :: Emission)
    sb = 5.6703744191844294e-8 # Stefan-Boltzmann constant
    
    θ = temperature
    θamb = emission.θamb
    θext = emission.θext

    h = emission.h
    ε₁ = emission.ε₁
    ε₂ = emission.ε₂

    @. flux = -h*(θ - θamb) - sb*(ε₁*θ^4 - ε₂*θext^4)
    return nothing
end






