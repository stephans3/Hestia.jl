#=
Author: Stephan Scholz
Year: 2021

This file contains a simulation of the heat equation
Dimension: 2D
Material: isotropic, dynamic
=#

using Hestia 
import LinearAlgebra

θ₀ = 600.0       # Initial temperature
λ = [10.0, 0.1]  # Thermal conductivity: temperature-dependend
ρ = [7800.0]     # Mass density: constant
c = [330.0, 0.4] # Specific heat capacity: temperature-dependend


L = 0.2    # Length
W = 0.2    # Width

Nx = 40    # Number of elements: x direction
Ny = 40    # Number of elements: y direction
Ntotal = Nx*Ny

property = createDynamicIsoProperty(λ, ρ, c)
plate    = HeatPlate(L, W, Nx, Ny, property)


### Boundaries ###
θamb = 300.0;

emission_nonlinear  = createEmission(10.0, 0.6, θamb)  # Nonlinear BC: heat transfer (linear) and heat radiation (quartic/nonlinear)

boundary_west   = :west 
boundary_east   = :east 
boundary_north  = :north  

boundary_plate = initBoundary(plate)
setEmission!(boundary_plate, emission_nonlinear, boundary_west)
setEmission!(boundary_plate, emission_nonlinear, boundary_east)
setEmission!(boundary_plate, emission_nonlinear, boundary_north)

const heatproblem = CubicHeatProblem(plate, boundary_plate)

function heat_conduction!(dθ, θ, param, t)
    property  = heatproblem.geometry.segmentation.heatProperty;
    boundary  = heatproblem.boundary

    diffusion!(dθ, θ, heatproblem.geometry, property, boundary)
end

# Initial conditions of ODE
θinit = θ₀*ones(Ntotal)
tspan = (0.0, 2000.0)
Δt    = 0.2             # Sampling time

import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)

using Plots
heatmap(reshape(sol[end], Nx, Ny))
