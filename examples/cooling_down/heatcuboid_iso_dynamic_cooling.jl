#=
Author: Stephan Scholz
Year: 2022

This file contains a simulation of the heat equation
Dimension: 3D
Material: isotropic, dynamic
=#

using Hestia 
import LinearAlgebra

θ₀ = 600.0  # Initial temperature
λ = [10.0, 0.1]  # Thermal conductivity: temperature-dependend
ρ = [7800.0]     # Mass density: constant
c = [330.0, 0.4] # Specific heat capacity: temperature-dependend

L = 0.2    # Length
W = 0.2    # Width    
H = 0.1    # Height

Nx = 24     # Number of elements: x direction
Ny = 24     # Number of elements: y direction
Nz = 10
Ntotal = Nx*Ny*Nz

property = createDynamicIsoProperty(λ, ρ, c)
cuboid   = HeatCuboid(L, W, H,Nx, Ny, Nz, property)


### Boundaries ###
θamb = 300.0;

emission_nonlinear = createEmission(10.0, 0.6, θamb)  # Nonlinear BC: heat transfer (linear) and heat radiation (quartic/nonlinear)
boundary_west  = :west 
boundary_east  = :east 
boundary_north = :north  
boundary_topside = :topside

boundary_cuboid = initBoundary(cuboid)
setEmission!(boundary_cuboid, emission_nonlinear, boundary_west)
setEmission!(boundary_cuboid, emission_nonlinear, boundary_east)
setEmission!(boundary_cuboid, emission_nonlinear, boundary_north)
setEmission!(boundary_cuboid, emission_nonlinear, boundary_topside)


function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, cuboid, property, boundary_cuboid)
end

# Initial conditions of ODE
θinit = θ₀*ones(Ntotal)
tspan = (0.0, 200.0)
Δt    = 0.2             # Sampling time

import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)

using Plots
heatmap(reshape(sol[end], Nx, Ny,Nz)[:,:,1])

heatmap(reshape(sol[end], Nx, Ny,Nz)[:,:,Nz])