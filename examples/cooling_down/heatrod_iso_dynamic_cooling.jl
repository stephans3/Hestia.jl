#=
Author: Stephan Scholz
Year: 2021

This file contains a simulation of the heat equation
Dimension: 1D
Material: isotropic, dynamic
=#

using Hestia 

θ₀ = 600.0      # Initial temperature
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

const heatproblem = CubicHeatProblem(heatrod, boundary_rod)

function heat_conduction!(dθ, θ, param, t)
    property  = heatproblem.geometry.segmentation.heatProperty;
    boundary  = heatproblem.boundary

    diffusion!(dθ, θ, heatproblem.geometry, property, boundary)
end

# Initial conditions of ODE
θinit = θ₀*ones(Nx)
tspan = (0.0, 2000.0)
Δt    = 0.2             # Sampling time

import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)

using Plots
heatmap(sol[:,:])