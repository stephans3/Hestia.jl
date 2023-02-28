#=
Author: Stephan Scholz
Year: 2022

This file contains a simulation of the heat equation with input
Dimension: 2D
Material: isotropic, dynamic
=#

using Hestia 

θ₀ = 300.0       # Initial temperature
λx = [10.0, 0.1]  # Thermal conductivity: temperature-dependend
λy = [5.0, 0.2]  # Thermal conductivity: temperature-dependend
ρ = [7800.0]     # Mass density: constant
c = [330.0, 0.4] # Specific heat capacity: temperature-dependend

L = 0.2    # Length
W = 0.2    # Width

Nx = 40    # Number of elements: x direction
Ny = 40    # Number of elements: y direction
Ntotal = Nx*Ny

property = createDynamicAnisoProperty(λx, λy, ρ, c)
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

### Actuation ###
num_actuators = 5        # Number of actuators per boundary
pos_actuators1 = :south  # Position of actuators
pos_actuators2 = :west   # Position of actuators

plate_actuation = initIOSetup(plate)

# Create actuator characterization
scale     = 1.0;
power     = 3;
curvature = 100.0;

config  = setConfiguration(scale, power, curvature)

setIOSetup!(plate_actuation, plate, num_actuators, config,  pos_actuators1)
setIOSetup!(plate_actuation, plate, num_actuators, config,  pos_actuators2, start_index = 6)


function heat_conduction!(dθ, θ, param, t)
    u_in = 4e5 * ones(2*num_actuators)    # heat input

    diffusion!(dθ, θ, plate, property, boundary_plate, plate_actuation, u_in)
end

# Initial conditions of ODE
θinit = θ₀*ones(Ntotal)
tspan = (0.0, 200.0)
Δt    = 0.2             # Sampling time

import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)

# Euler method
# sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)

# Runge-Kutta method
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.KenCarp5(), saveat=1.0)

using Plots
heatmap(reshape(sol[end], Nx, Ny))