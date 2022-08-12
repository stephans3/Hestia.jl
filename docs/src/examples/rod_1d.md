# One-dimensional rod

In this example we assume *isotropic*, temperature-*independent* heat conduction in a 1-dimensional rod. The rod is heated-up on the left side ($x=0$) and emits thermal energy on the right side ($x=L$) via convection.  The heat equation is noted as 

$\frac{\partial \theta(t,x)}{\partial t} = \frac{\lambda}{c~\rho} \frac{\partial^2 \theta(t,x)}{\partial x^2}$

with boundary conditions

$-\lambda \frac{\partial \theta(t,x)}{\partial x} = b ~ u(t)$

which is the induced heat flux on the left side $x=0$ and

$\lambda \frac{\partial \theta(t,x)}{\partial x} = -h (\theta(t,x) - \theta_{amb})$

which is the emitted heat flux on the right side $x=L$. We assume an initial and ambient temperature of $300$ Kelvin. 

All used values are listed in the table below.

| Name                      | Variable      | Value / Formula     |
| :---                      |    :----:     |  :--- |
| Length                    | L             | 0.2   |
| Discretization            | $N_{x}$       | 40    |
| Density                   | $\rho$        | 8000  |
| Thermal conductivity      | $\lambda$     | 50    |
| Specific heat capacity    | $c$           | 400   |
| Heat transfer coefficient | $h$           | 5     |
| Ambient temperature       | $\theta_{amb}$| 300   |
| Spatial characterization  | $b$           | 1     |


# 1. Material Properties and Geometry
In the first step we define the material properties and the geometry. As stated in the beginning we assume *isotropic*, temperature-*independent* material properties and thus we create a `StaticIsoProperty` with `createStaticIsoProperty(λ, ρ, c)`. Next, we create the geometry: here the one-dimensional rod `HeatRod`.

```julia
λ = 50.0    # Thermal conductivity
ρ = 8000.0  # Mass density
c = 400.0   # Specific heat capacity
property = createStaticIsoProperty(λ, ρ, c)

L = 0.2     # Length
Nx = 40     # Number of elements: x direction
heatrod  = HeatRod(L, Nx, property)
```

# 2. Emission
Now, we define the emitted heat flux on the right boundary side ($x=L$). On the right boundary side we only assume linear heat transfer (convection) and no heat radiation, and so emissivity $\epsilon=0$. 

```julia
h = 5.0       # Heat transfer coefficient
ϵ = 0.0       # Emissivity
θamb = 300.0; # Ambient temperature
emission  = createEmission(h, ϵ, θamb)   # Only convection
```
Next, this emission has to be assigned for the right boundary side, which is specified as `:east`. To do so, we initialize the boundary with `initBoundary()` to yield a container which stores all emissions for each boundary.

```julia
boundary_rod  = initBoundary(heatrod)
boundary_side = :east 
setEmission!(boundary_rod, emission, boundary_side)
```

# 3. Actuation
Before we define the actuation we initialize the `IOSetup` which is a container for actuator and sensor setups.

```julia
rod_actuation = initIOSetup(heatrod)
```
In this example, we assume one input signal on the left boundary side (`:west`). In this example, the spatial characterization is only a single point described by the scaling value. Finally, the actuator setup on the left boundary is created with `setIOSetup!()`. 

```julia
scale         = 1.0;  # b=1
input_id      = 1     # Actuator index
pos_actuators = :west # Position of actuators
setIOSetup!(rod_actuation, heatrod, input_id, scale,  pos_actuators)
```

# 4. Simulation

The simulation is built as explained in the [DifferentialEquations.jl documentation](https://diffeq.sciml.ai/stable/). We define the differential equation interface `heat_conduction!(dθ, θ, param, t)` with temperature `θ` and use method `diffusion!` to simulate the diffusion. We assume a constant input $u(t) = 4 \cdot 10^5$. Only vector-valued input signals are supported and so it is multiplied with `ones(1)`.

```julia
function heat_conduction!(dθ, θ, param, t)
    u_in = 4e5 * ones(1)    # heat input as vector
    diffusion!(dθ, θ, heatrod, property, boundary_rod, rod_actuation, u_in)
end
```

As noted in the beginning we assume an initial temperature of $300$ Kelvin for the whole rod, and we simulate the diffusion for $600$ seconds.

```julia
θinit = 300*ones(Nx)
tspan = (0.0, 600.0)
```
### Forward Euler method
If we use the forward Euler method we have to fix the sampling time `Δt` such that the numerical integration is stable check the [numerical stability](https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis). From **von Neumann stability analysis** we know that the condition

$\Delta t < \frac{\Delta x^2}{2} \left(\frac{\lambda}{c~\rho}\right)^{-1}$ 

has to hold to guarantee numerical stability. 

```julia
Δx = heatrod.sampling[1]  # Spatial discretization
α = λ/(ρ *c)              # Diffusivity constant
dt_max = Δx^2/(2*α)       # Highest possible sampling time
```
Here, the sampling time has to be greater than $0.8$ seconds and thus we choose the sampling time `Δt = 0.5`. In the final step, we define numerical integration method (here: `Euler()`), the `ODEProblem` and solve the differential equation. Here, we save the integration for every second `saveat=1.0`.

```julia
import OrdinaryDiffEq
alg = OrdinaryDiffEq.Euler();
prob = OrdinaryDiffEq.ODEProblem(heat_conduction!,θinit,tspan)
sol = OrdinaryDiffEq.solve(prob,alg, dt=Δt, saveat=1.0)
```

### Adaptive numerical integration
It is also possible to use adaptive integration methods like `KenCarp5()`. Then, we do not set a fixed sampling time and we do not have to check the numerical stability.

```julia
alg_2 = OrdinaryDiffEq.KenCarp5()
sol_2 = OrdinaryDiffEq.solve(prob,alg_2, saveat=1.0)
```

## Additional: Proportional control
In this section we want to implement a proportional controller for the one-dimensional heat equation. In general proportional controllers are defined by

$u(t) = K_{p} ~ (r(t) - y(t))$

with proportional gain $K_{p}$, reference $r(t)$ and system output $y(t)$. Here, we specify the right boundary as the system output $y(t) = \theta(t,L)$  (`θ[end]`) and we assume a constant reference value $r=400$ Kelvin. The proportional gain is found manually as $K_{p} = 600$. In this example we assume that our actuator is only able to heat - not to cool. This means, the control signal has to be positive or zero and thus we change the control law to

$u(t) = K_{p} ~ \max(r(t) - y(t), 0).$

```julia
Kp = 600    # Proportional gain
Θref = 400  # Reference temperature of right boundary
controller(ref, yout) = Kp*max((ref-yout),0) # Proportional controller
```

In the differential equations interface we use the `controller()` method to calculate the input signal.

```julia
function heat_conduction_controlled!(dθ, θ, param, t)
    u_in = controller(Θref, θ[end]) * ones(1)    # heat input
    diffusion!(dθ, θ, heatrod, property, boundary_rod, rod_actuation, u_in)
end
```

Finally, we use a longer time span up to $3500$ seconds because the controller is not able to reach the reference in $600$ seconds as above.  

```julia
tspan_cntr = (0.0, 3500.0)
prob_cntr = OrdinaryDiffEq.ODEProblem(heat_conduction_controlled!,θinit,tspan_cntr)
sol_cntr = OrdinaryDiffEq.solve(prob_cntr,alg, dt=Δt, saveat=1.0)
```

# Full source code listing


```julia
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


### Forward Euler method ###
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


### Adaptive numerical integration ###
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


### Create plots ###
using Plots
heatmap(sol[:,:], title="1D Heat conduction with Euler method")
heatmap(sol_2[:,:], title="1D Heat conduction with KenCarp5")


### Plots of proportional control ###
plot(sol_cntr,legend=false, title="Proportional control of 1D heat conduction")

u_signals = controller.(Θref, sol_cntr[end,1:end-1])
tgrid = 0 : 1.0 : tspan_cntr[2]-1
plot(tgrid, u_signals, title="Control signals of proportional control")
```

