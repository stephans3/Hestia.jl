#=
Author: Stephan Scholz
Year: 2021

This file contains a simulation of the heat equation
Dimension: 1D
Material: isotropic, static
=#

using Hestia 

θ₀ = 600.0  # Initial temperature
λ = 45.0    # Thermal conductivity: constant
ρ = 7800.0  # Mass density: constant
c = 480.0   # Specific heat capacity: constant

L = 0.2     # Length
Nx = 40    # Number of elements: x direction

property = createStaticIsoProperty(λ, ρ, c)
heatrod  = HeatRod(L, Nx, property)

### Boundaries ###
θamb = 300.0;

emission_right  = createEmission(5.0, 0.5, θamb)   # Nonlinear BC: heat transfer (linear) and heat radiation (quartic/nonlinear)
boundary_right  = :east 

boundary_rod = initBoundary(heatrod)
setEmission!(boundary_rod, emission_right, boundary_right)

function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, heatrod, property, boundary_rod)
end

# Initial conditions of ODE
θinit = θ₀*ones(Nx)
tspan = (0.0, 2000.0)
Δt    = 0.2             # Sampling time

Δx = heatrod.sampling[1]
α = λ/(ρ *c)
min_dt = 0.5*Δx^2/α

if Δt > min_dt
    error("Numerical stability is not guaranteed! Choose a smaller sampling time.")
end

import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)

using Plots
heatmap(sol[:,:])
