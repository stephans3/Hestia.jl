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

L = 0.2    # Length
W = 0.2    # Width    
H = 0.1    # Height

Nx = 24     # Number of elements: x direction
Ny = 24     # Number of elements: y direction
Nz = 10
Ntotal = Nx*Ny*Nz

property = createStaticIsoProperty(λ, ρ, c)
cuboid   = HeatCuboid(L, W, H, Nx, Ny, Nz, property)

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

### Actuation ###
num_actuators = (4,3)        # Number of actuators
pos_actuators = :underside   # Position of actuators
num_act_total = num_actuators[1]*num_actuators[2]

# cuboid_actuation = initActuation(cuboid)
cuboid_actuation = initIOSetup(cuboid)

# Create actuator characterization
scale     = 1.0;
power     = 3;
curvature = 100.0;

config  = setConfiguration(scale, power, curvature)

# setActuation!(cuboid_actuation, cuboid, num_actuators, config,  pos_actuators)
setIOSetup!(cuboid_actuation, cuboid, num_actuators, config,  pos_actuators)


function heat_conduction!(dθ, θ, param, t)
    u_in = 4e5 * ones(num_act_total)

    diffusion!(dθ, θ, cuboid, property, boundary_cuboid, cuboid_actuation, u_in)
end

# Initial conditions of ODE
θinit = θ₀*ones(Ntotal)
tspan = (0.0, 200.0)
Δt    = 0.2             # Sampling time

# Check numerical stability for Euler forward method
Δx = cuboid.sampling[1]
Δy = cuboid.sampling[2]
Δz = cuboid.sampling[3]
α = λ/(ρ *c)
min_dt = 0.5*inv(1/Δx^2 + 1/Δy^2 + 1/Δz^2)/α

if Δt > min_dt
    error("Numerical stability is not guaranteed! Choose a smaller sampling time.")
end

import LinearAlgebra
import OrdinaryDiffEq

prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
# sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.KenCarp5(), saveat=1.0)

using Plots
heatmap(reshape(sol[end], Nx, Ny,Nz)[:,:,1])

heatmap(reshape(sol[end], Nx, Ny,Nz)[:,:,Nz])
