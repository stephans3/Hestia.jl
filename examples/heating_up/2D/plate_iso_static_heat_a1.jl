#=
Author: Stephan Scholz
Year: 2022

This file contains a simulation of the heat equation with input
Dimension: 2D
Material: isotropic, static
=#

using Hestia 

θ₀ = 300.0  # Initial temperature
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

emission_robin    = createEmission(10.0, 0.6, θamb)  # Nonlinear BC: heat transfer (linear) and heat radiation (quartic/nonlinear)

boundary_west = :west 
boundary_east = :east 
boundary_north = :north  

boundary_plate = initBoundary(plate)
setEmission!(boundary_plate, emission_robin, boundary_west)
setEmission!(boundary_plate, emission_robin, boundary_east)
setEmission!(boundary_plate, emission_robin, boundary_north)

### Actuation ###
num_actuators = 5        # Number of actuators
pos_actuators = :south   # Position of actuators

# plate_actuation = initActuation(plate)
plate_actuation = initIOSetup(plate)


# Create actuator characterization
scale     = 1.0;
power     = 3;
curvature = 100.0;

config  = setConfiguration(scale, power, curvature)

# setActuation!(plate_actuation, plate, num_actuators, config,  pos_actuators)
setIOSetup!(plate_actuation, plate, num_actuators, config,  pos_actuators)


function heat_conduction!(dθ, θ, param, t)
    u_in = 4e5 * ones(num_actuators)    # heat input

    diffusion!(dθ, θ, plate, property, boundary_plate, plate_actuation, u_in)
end

θinit = θ₀*ones(Ntotal)

tspan = (0.0, 200.0)
Δt = 1e-2               # Sampling time


import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
# sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.KenCarp5(), saveat=1.0)

using Plots
heatmap(reshape(sol[end], Nx, Ny))
