using Hestia
λ = [20.0, 0.1]  # Thermal conductivity: temperature-dependend
ρ = [8000.0]     # Mass density: constant
c = [220.0, 0.6] # Specific heat capacity: temperature-dependend
property = DynamicIsotropic(λ, ρ, c)

# Geometry
L = 0.2    # Length
W = 0.1    # Width
Nx = 60    # Number of elements: x direction
Ny = 40    # Number of elements: y direction
plate = HeatPlate(L, W, Nx, Ny)

# Emission
h = 5.0       # Heat transfer coefficient
ϵ = 0.3       # Emissivity
θamb = 300.0; # Ambient temperature
emission = Emission(h, θamb, ϵ)   # Convection and Radiation

# Boundary Conditions
boundary  = Boundary(plate)
setEmission!(boundary, emission, :east)
setEmission!(boundary, emission, :west)
setEmission!(boundary, emission, :north)

# Actuation
actuation = IOSetup(plate)
pos_actuators = :south  # Position of actuators

# Actuator b1 and b2
config1  = RadialCharacteristics(1, 1, 4) # Set actuators characterization b1, b2
config2  = RadialCharacteristics(0.9, 2, 20) # Set actuators characterization b3

partition    = [1,2,3]
config_table = [config1, config1, config2]
setIOSetup!(actuation, plate, partition, config_table,  pos_actuators)


### Simulation ###
function heat_conduction!(dθ, θ, param, t)
    u_in = 2e5 * ones(3)    # heat input
    diffusion!(dθ, θ, plate, property, boundary, actuation, u_in)
end

θinit = 300*ones(Nx*Ny) # Inital temperatures
tspan = (0.0, 300.0) # Time span

using OrdinaryDiffEq
alg = KenCarp5()
prob = ODEProblem(heat_conduction!,θinit,tspan)
sol = solve(prob,alg, saveat=1.0)

xgrid = L/(2Nx) : L/Nx : L # Position in x-direction
ygrid = W/(2Ny) : W/Ny : W # Position in y-direction

using Plots
act_char, _ = getCharacteristics(actuation, pos_actuators)
spat_char = scatter(xgrid, act_char, xlabel="Length L", ylabel="Spatial Characteristics",label=["Actuator 1" "Actuator 2" "Actuator 3"], legend=:bottomright, linewidth=3)

anim = @animate for i in axes(sol,2)
    heatmap(xgrid, ygrid, reshape(sol[:,i], Nx, Ny)', xlabel="Length L", ylabel="Width W", title=string("Temperature at Time ",sol.t[i], " seconds") )
end
gif(anim, "plate_2d_animation.gif", fps=10)

savefig(spat_char, "docs/assets/spatial_characteristics.png")