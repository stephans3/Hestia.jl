#=
Author: Stephan Scholz
Year: 2022

This file contains a simulation of the heat equation with input
Dimension: 2D
Material: isotropic, static
=#

using Hestia 

θ₀ = 300.0   # Initial temperature
λ  = 45.0    # Thermal conductivity: constant
ρ  = 7800.0  # Mass density: constant
c  = 480.0   # Specific heat capacity: constant

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

# const heatproblem = CubicHeatProblem(plate, boundary_plate)


### Actuation ###
num_actuators1  = 5        # Number of actuators
num_actuators2  = 7        # Number of actuators
pos_actuators1 = :south   # Position of actuators
pos_actuators2 = :north    # Position of actuators

# plate_actuation = initActuation(plate)
plate_actuation = initIOSetup(plate)

# Create actuator characterization
scale     = 1.0;
power     = 3;
curvature = 100.0;

config  = setConfiguration(scale, power, curvature)
config_table1 = RadialConfiguration[]
config_table2 = RadialConfiguration[]

for i = 1 : num_actuators1
    config_table1 = vcat(config_table1, config)
end

for i = 1 : num_actuators2
    config_table2 = vcat(config_table2, config)
end

#=
Partitions are used to define which input signals act at which position.
A partition entry zero means that this partition cell corresponds to no input signal.
One input signal can be used for more than one partition cell if partition entry has the same value.
=#
partition1 = collect(1:num_actuators1)
partition1[2] = 0

partition2 = collect(num_actuators1+1:num_actuators1+num_actuators2)
partition2[3] = 0
partition2[5] = 0


# setActuation!(plate_actuation, plate, partition1, config_table1,  pos_actuators1)
# setActuation!(plate_actuation, plate, partition2, config_table2,  pos_actuators2)

setIOSetup!(plate_actuation, plate, partition1, config_table1,  pos_actuators1)
setIOSetup!(plate_actuation, plate, partition2, config_table2,  pos_actuators2)


# Uncontrolled / free system
function heat_conduction!(dθ, θ, param, t)
    # property = heatproblem.geometry.segmentation.heatProperty;
    # boundary  = heatproblem.boundary

    u_in = zeros(num_actuators1 + num_actuators2)
    u_in[1:num_actuators1] .= 4e5
    u_in[num_actuators1+1:end] .= 3e5

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
