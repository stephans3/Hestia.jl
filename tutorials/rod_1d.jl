using Hestia 

λ = 50.0    # Thermal conductivity: constant
ρ = 8000.0  # Mass density: constant
c = 400.0   # Specific heat capacity: constant
property = createStaticIsoProperty(λ, ρ, c)

L = 0.2     # Length
Nx = 40    # Number of elements: x direction
heatrod  = HeatRod(L, Nx, property)

### Boundaries ###
h = 5.0
ϵ = 0.0
θamb = 300.0;
emission  = createEmission(h, ϵ, θamb)   # Only convection

boundary_rod  = initBoundary(heatrod)
boundary_side = :east 
setEmission!(boundary_rod, emission, boundary_side)


### Actuation ###
rod_actuation = initIOSetup(heatrod)

scale         = 1.0;  # b=1
input_id      = 1     # Actuator index
pos_actuators = :west # Position of actuators
setIOSetup!(rod_actuation, heatrod, input_id, scale,  pos_actuators)



### Simulation ###
function heat_conduction!(dθ, θ, param, t)
    u_in = 4e5 * ones(1)    # heat input
    diffusion!(dθ, θ, heatrod, property, boundary_rod, rod_actuation, u_in)
end

θinit = 300*ones(Nx)
tspan = (0.0, 600.0)

### Forward Euler method
Δx = heatrod.sampling[1]  # Spatial discretization
α = λ/(ρ *c)              # Diffusivity constant
dt_max = Δx^2/(2*α)       # Highest possible sampling time

Δt = 0.5

if Δt > dt_max
    error("Numerical stability is not guaranteed! Choose a smaller sampling time.")
end

import OrdinaryDiffEq
alg = OrdinaryDiffEq.Euler();
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,alg, dt=Δt, saveat=1.0)


### Adaptive numerical integration
alg_2 = OrdinaryDiffEq.KenCarp5()
sol_2 = OrdinaryDiffEq.solve(prob,alg_2, saveat=1.0)


### Proportional control ###
Kp = 600    # Proportional gain
Θref = 400  # Reference temperature of right boundary
controller(ref, yout) = Kp*max((ref-yout),0) # Proportional controller

function heat_conduction_controlled!(dθ, θ, param, t)
    u_in = controller(Θref, θ[end]) * ones(1)    # heat input
    diffusion!(dθ, θ, heatrod, property, boundary_rod, rod_actuation, u_in)
end

tspan_cntr = (0.0, 3500.0)
prob_cntr = OrdinaryDiffEq.ODEProblem(heat_conduction_controlled!,θinit,tspan_cntr)
sol_cntr = OrdinaryDiffEq.solve(prob_cntr,alg, dt=Δt, saveat=1.0)


### Create plots
using Plots
heatmap(sol[:,:], title="1D Heat conduction with Euler method")
heatmap(sol_2[:,:], title="1D Heat conduction with KenCarp5")


## Plots of proportional control
plot(sol_cntr,legend=false, title="Proportional control of 1D heat conduction")

u_signals = controller.(Θref, sol_cntr[end,1:end-1])
tgrid = 0 : 1.0 : tspan_cntr[2]-1
plot(tgrid, u_signals, title="Control signals of proportional control")
