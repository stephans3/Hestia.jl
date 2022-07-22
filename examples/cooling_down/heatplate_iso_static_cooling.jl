#=
Author: Stephan Scholz
Year: 2021

This file contains a simulation of the heat equation
Dimension: 2D
Material: isotropic, static
=#

using Hestia 

θ₀ = 600.0  # Initial temperature for all cells
λ = 45.0    # Thermal conductivity: constant
ρ = 7800.0  # Mass density: constant
c = 480.0   # Specific heat capacity: constant

L = 0.2     # Length
W = 0.2     # Width

Nx = 40     # Number of elements: x direction
Ny = 40     # Number of elements: y direction
Ntotal = Nx*Ny

property = createStaticIsoProperty(λ, ρ, c)
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

function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, plate, property, boundary_plate)
end

# Initial conditions of ODE
θinit = θ₀*ones(Ntotal)
tspan = (0.0, 2000.0)
Δt    = 0.2             # Sampling time

# Check numerical stability for Euler forward method
Δx = plate.sampling[1]
Δy = plate.sampling[2]
α = λ/(ρ *c)
min_dt = 0.5*inv(1/Δx^2 + 1/Δy^2)/α

if Δt > min_dt
    error("Numerical stability is not guaranteed! Choose a smaller sampling time.")
end

import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)


using Plots
heatmap(reshape(sol[end], Nx, Ny))
