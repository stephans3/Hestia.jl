#=
Author: Stephan Scholz
Year: 2022

This file contains a simulation of the heat equation with input
Dimension: 1D
Material: isotropic, dynamic
=#

using Hestia 

θ₀ = 300.0  # Initial temperature
λ = [8.0, 0.1]  # Thermal conductivity
ρ = [7800.0]    # Mass density: constant
c = [330, 0.5]  # Specific heat capacity

L = 0.2     # Length
Nx = 40    # Number of elements: x direction

property = createDynamicIsoProperty(λ, ρ, c)
heatrod  = HeatRod(L, Nx, property)


### Boundaries ###
θamb = 300.0;

emission_right  = createEmission(5.0, 0.5, θamb)   # Nonlinear BC: heat transfer (linear) and heat radiation (quartic/nonlinear)
boundary_right  = :east 

boundary_rod = initBoundary(heatrod)
setEmission!(boundary_rod, emission_right, boundary_right)


### Actuation ###
pos_actuators = :west   # Position of actuators

rod_actuation = initIOSetup(heatrod)

# Create actuator characterization
scale     = 0.9;
input_id = 1

setIOSetup!(rod_actuation, heatrod, input_id, scale,  pos_actuators)

function heat_conduction!(dθ, θ, param, t)
    u_in = 4e5 * ones(1)    # heat input

    diffusion!(dθ, θ, heatrod, property, boundary_rod, rod_actuation, u_in)
end

# Initial conditions of ODE
θinit = θ₀*ones(Nx)
tspan = (0.0, 200.0)
Δt    = 0.2             # Sampling time

import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)

# Euler method 
# sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)

# Runge-Kutta method
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.KenCarp5(), saveat=1.0)


using Plots
heatmap(sol[:,:])
