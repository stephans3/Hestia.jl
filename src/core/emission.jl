abstract type AbstractEmission end


#############################################################

# Constant heat radiation + conduction

#############################################################


"""
    Emission  <: AbstractEmission

Returns the right-hand side of a natural Robin boundary condition including heat conduction and emission to the environment.

### Elements

`conduction` : heat conduction coefficient

`radiation` : heat radiation coefficient

`θamb` : ambient temperature
"""
mutable struct Emission <: AbstractEmission 
   conduction :: Real  # Heat conduction coefficient
   radiation :: Real   # Heat radiation coefficient
   θamb :: Real          # Ambient temperature  
end



"""
    createEmission(conduction :: Real, emissivity :: Real, θamb :: Real)

Returns an `Emission` for a given `conductivity`, `emissivity` and ambient temperature `θamb`

### Information
This constant is used:

- Stefan-Boltzmann constant: `σ = 5.6703744191844294e-8`

"""
function createEmission(conduction :: Real, emissivity :: Real, θamb :: Real)
    if emissivity < 0 || emissivity > 1
        error("Emissivity has to be in interval [0, 1]!")
    elseif conduction < 0
        error("Conduction has to be greater than zero!")
    end

    
    sb = 5.6703744191844294e-8

    radiation = emissivity * sb

    return Emission(conduction, radiation, θamb)
end




"""
    emit!(flux :: Array{S1,1} where S1 <: Real, temperature :: Array{S2,1} where S2 <: Real, emission :: Emission)

Calculates the right-hand side of the natural Robin boundary along a boundary for a given `Emission`.

Note: in-place operation - results are saved in array `flux`.
"""
function emit!(flux :: Array{S1,1} where S1 <: Real, temperature :: Array{S2,1} where S2 <: Real, emission :: Emission)
    θ = temperature
    θamb = emission.θamb
    h = emission.conduction
    k = emission.radiation

    flux .= (  -h *( θ .- θamb) - k*( θ.^4 .- θamb^4)  )

    return nothing
end


"""
    emit!(flux :: Array{S1,1} where S1 <: Real, temperature :: Array{S2,1} where S2 <: Real, emission :: Emission)

Calculates the right-hand side of the natural Robin boundary along a boundary for a given `Emission`.

Note: in-place operation - results are saved in array `flux`.
"""
function emit(temperature :: Real, emission :: Emission)
    θ = temperature
    θamb = emission.θamb
    h = emission.conduction
    k = emission.radiation

    flux = (  -h *( θ - θamb) - k*( θ^4 - θamb^4)  )

    return flux
end



