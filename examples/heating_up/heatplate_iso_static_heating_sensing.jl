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

const heatproblem = CubicHeatProblem(plate, boundary_plate)


### Actuation ###
num_actuators = 5        # Number of actuators
pos_actuators = :south   # Position of actuators

plate_actuation = initIOSetup(plate)

# Create actuator characterization
scale     = 1.0;
power     = 3;
curvature = 100.0;

config  = setConfiguration(scale, power, curvature)

setIOSetup!(plate_actuation, plate, num_actuators, config,  pos_actuators)

function heat_conduction!(dθ, θ, param, t)
    property  = heatproblem.geometry.segmentation.heatProperty;
    boundary  = heatproblem.boundary

    u_in = 4e7 * ones(num_actuators) # 4e5 * ones(num_actuators)    # heat input

    diffusion!(dθ, θ, heatproblem.geometry, property, boundary, plate_actuation, u_in)
    # diffusion!(dθ, θ, heatproblem.geometry, property, boundary)
end

θinit = θ₀*ones(Ntotal)# θ₀*(1 .+ 0.1*sin.(1:1:Ntotal)) # θ₀*ones(Ntotal)

tspan = (0.0, 200.0)
Δt = 1e-2               # Sampling time



import OrdinaryDiffEq
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.KenCarp5(), dt=Δt, saveat=1.0)
# sol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Euler(), dt=Δt, saveat=1.0)


### Sensor ###
num_sensor = 3        # Number of sensors
pos_sensor = :north   # Position of sensors

plate_sensing = initIOSetup(plate)

# Create actuator characterization
scale     = 1.0;
power     = 2;
curvature = 100.0;

config_sensor  = setConfiguration(scale, power, curvature)

setIOSetup!(plate_sensing, plate, num_sensor, config,  pos_sensor)

measurementData = zeros(201,0)

for i = 1 : num_sensor
    cellidx, chars = getSensing(plate_sensing, i, pos_sensor)
    temps = sol[cellidx,1:end]
    measure = measureWAM(temps, chars)
    measurementData = hcat(measurementData, measure)
end

using Plots
plot(sol.t, measurementData, title="Measured output", xlabel="Time [s]", ylabel="Temperature [K]", legend=:bottomright)

