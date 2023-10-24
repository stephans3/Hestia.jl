abstract type AbstractEmission end


#############################################################

# Constant heat radiation + conduction

#############################################################


"""
    Emission  <: AbstractEmission

Type Emission contains the coefficients for heat transfer (convection) `h`, heat radiation `k`, and the ambient temperature `θamb` to model linear or nonlinear Stefan-Boltzmann boundary conditions
    
```math
\\Phi = -h ~ (\\theta - \\theta_{amb}) -k ~ (\\theta^4 - \\theta_{amb}^4)
```    

Constructor `Emission(h,ϵ,θamb)` expects emissivity `ϵ` which must be in interval [0,1]. The heat radiation coefficient is calculated internally as `k=ϵ⋅σ` using the Stefan-Boltzmann constant: `σ = 5.6703744191844294e-8`.

### Elements

`h` : heat transfer coefficient

`k` : heat radiation coefficient

`θamb` : ambient temperature
"""
mutable struct Emission <: AbstractEmission 
   h :: Real  # Heat conduction coefficient
   k :: Real   # Heat radiation coefficient
   θamb :: Real        # Ambient temperature  

   function Emission(h, ϵ, θamb)
    if ϵ < 0 || ϵ > 1
        error("Emissivity has to be in interval [0, 1]!")
    elseif h < 0
        error("Heat transfer coefficient has to be greater than zero!")
    end
    sb = 5.6703744191844294e-8
    new(h, sb*ϵ, θamb )
   end
end




"""
    emit(temperature :: Real, emission :: Emission)

Calculates the right-hand side of the boundary conditions for a given `Emission`.

```math
\\Phi = -h ~ (\\theta - \\theta_{amb}) -k ~ (\\theta^4 - \\theta_{amb}^4)
```
"""
function emit(temperature :: Real, emission :: Emission)
    θ = temperature
    θamb = emission.θamb
    h = emission.h
    k = emission.k

    flux = (  -h *( θ - θamb) - k*( θ^4 - θamb^4)  )

    return flux
end


"""
    emit!(flux :: Vector{<:Real}, temperature :: Vector{<:Real}, emission :: Emission)

Calculates the right-hand side of the natural Robin boundary along a boundary for a given `Emission`.

Note: in-place operation - results are saved in array `flux`.
"""
function emit!(flux :: Vector{<:Real}, temperature :: Vector{<:Real}, emission :: Emission)
    θ = temperature
    θamb = emission.θamb
    h = emission.h
    k = emission.k

    @. flux = (-h*(θ - θamb) - k*(θ^4 - θamb^4))
    return nothing
end






