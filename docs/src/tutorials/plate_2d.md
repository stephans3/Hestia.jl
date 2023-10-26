# Two-dimensional plate
In this example we assume isotropic, temperature-*dependend* heat conduction in a 2-dimensional plate ($(x,y) \in [0,L] \times [0,W]$) with length $L=0.2$ and width $W=0.1$. The plate is heated-up with 3 actuators on boundary 
- `:south` ($(x,y) \in [0,L] \times \{0\}$) 

and emits thermal energy on boundaries 
- `:west` - $(x,y) \in \{0\} \times [0,W]$, 
- `:east` - $(x,y) \in \{L\} \times [0,W]$ and 
- `:north` - $(x,y) \in [0,L] \times \{W\}$.  

The quasi-linear heat equation is noted as

$\frac{\partial \theta(t,x,y)}{\partial t} = \frac{1}{c(\theta)~\rho(\theta)} \left[\frac{\partial}{\partial x} \left( \lambda(\theta) \frac{\partial \theta(t,x)}{\partial x}\right) + \frac{\partial}{\partial y} \left( \lambda(\theta) \frac{\partial \theta(t,x,y)}{\partial y}\right) \right]$


with boundary conditions
$-\lambda \frac{\partial \theta(t,x,y)}{\partial y} = b_{1}(x,y) ~ u_{1}(t) + b_{2}(x,y) ~ u_{2}(t) + b_{3}(x,y) ~ u_{3}(t)$

as actuation on boundary `:south` at $y=0$ and emissions

$-\lambda \frac{\partial \theta(t,x,y)}{\partial x} = -h (\theta(t,x,y) - \theta_{amb}) - \epsilon \sigma (\theta(t,x,y)^4 - \theta_{amb}^4)$

on boundary `:west` at $x=0$

$\lambda \frac{\partial \theta(t,x,y)}{\partial x} = -h (\theta(t,x,y) - \theta_{amb}) - \epsilon \sigma (\theta(t,x,y)^4 - \theta_{amb}^4)$

on boundary `:east` at $x=L$ and 

$\lambda \frac{\partial \theta(t,x,y)}{\partial y} = -h (\theta(t,x,y) - \theta_{amb}) - k (\theta(t,x,y)^4 - \theta_{amb}^4)$

on boundary `:north` at $y=W$. The ambient temperature $\theta_{amb}$, the heat transfer coefficient $h$ and0 the heat radiation coefficient $k = \epsilon ~ \sigma$ with emissivity $\epsilon$  and Stefan-Boltzmann constant $\sigma \approx 5.67 \cdot 10^{-8} \frac{W}{m^2 K^4}$ are assumed to be equal for each boundary side. 

We assume an initial and ambient temperature of $300$ Kelvin. 

All used values are listed in the table below.

| Name                      | Variable      | Value / Formula     |
| :---                      |    :----:     |  :--- |
| Length                    | L             | 0.2   |
| Width                     | W             | 0.1   |
| Discretization            | $N_{x}$       | 40    |
|                           | $N_{y}$       | 20    |
| Density                   | $\rho$        | 8000              |
| Thermal conductivity      | $\lambda$     |  $20 + 0.1 \theta$  |
| Specific heat capacity    | $c$           | $220 + 0.6 \theta$  |
| Heat transfer coefficient | $h$           | 5     |
| Emissivity                | $\epsilon$    | 0.3   |
| Ambient temperature       | $\theta_{amb}$| 300   |

## Material Properties and Geometry
In the first step we define the material properties and the geometry. As stated in the beginning we assume *isotropic*, temperature-*dependent* thermal conductivity
```math
\lambda(\theta) = \lambda_{0} + \lambda_{1} \theta =  20 + 0.1 \theta
```
and the temperature-*dependent* specific heat capacity
```math
c(\theta) = c_{0} + c_{1} \theta =  220 + 0.6 \theta
```
and save the parameters in `DynamicIsotropic`. 

```@example 2d_plate
using Hestia
λ = [20.0, 0.1]  # Thermal conductivity: temperature-dependend
ρ = [8000.0]     # Mass density: constant
c = [220.0, 0.6] # Specific heat capacity: temperature-dependend
property = DynamicIsotropic(λ, ρ, c)
```

We create the geometry: here the two-dimensional plate `HeatPlate`.
```@example 2d_plate
L = 0.2    # Length
W = 0.1    # Width
Nx = 40    # Number of elements: x direction
Ny = 40    # Number of elements: y direction
plate    = HeatPlate(L, W, Nx, Ny)
```

## Emission
Next, we define the emitted heat flux on the boundary sides `:west`, `:east` and `:north`. We assume the same linear heat transfer (convection) and nonlinear heat radiation for all boundaries as noted in the table above. The constructor of `Emission` expects the heat transfer coefficient and the emissivity; the heat radiation coefficient is computed internally using the Stefan-Boltzmann constant. 

```@example 2d_plate
h = 5.0       # Heat transfer coefficient
ϵ = 0.3       # Emissivity
θamb = 300.0; # Ambient temperature
emission  = Emission(h, ϵ, θamb)   # Convection and Radiation
```
The emission has to be assigned for all boundary sides. We initialize the boundary with `Boundary()` to yield a container which stores all emissions and set the emissions. The emission at boundary `:south` is initiallized automatically as zero-Neumann boundary condition (thermal insulation) because only heat supply is modelled on this boundary side.

```@example 2d_plate
boundary  = Boundary(plate)
setEmission!(boundary, emission, :east)
setEmission!(boundary, emission, :west)
setEmission!(boundary, emission, :north)
```

# 3. Actuation
Firstly, we initialize `IOSetup` as a container for actuator and sensor setups. In this example, we assume three input signals on boundary side `:south`. 
```@example 2d_plate
actuation = IOSetup(plate)
pos_actuators = :south;  # Position of actuators
```

The spatial characteristics is defined by

$b(x) = m ~ \exp(-M ~ \lVert x - x_{c} \rVert^{2 \nu})$

with scaling $m \in [0,1]$, curvature matrix $M \in \mathbb{R}^{3 \times 3}$, power $\nu \in \{1,2,\cdots\}$ and central point $x_{c}$. Curvature matrix is a diagonal matrix and usually defined as a scaled identity. Central point $x_{c}$ is computed internally. The first both actuators (from left to right) have the same characteristics $(m_{1}, \nu_{1}, M_{1})$ and the third actuator has characteristics $(m_{2}, \nu_{2}, M_{2})$ as listed in the Table below.


| Name               | Variable                 | Value / Formula              |
| :---               |    :----:                |  :---                        |
| Actuator b1, b2    | $(m_{1}, \nu_{1}, M_{1})$|  $(1.0, 1, 4~I_{3\times3})$  |
| Actuator b3        | $(m_{2}, \nu_{2}, M_{2})$|  $(0.9, 2, 20~I_{3\times3})$ |

We define two characteristics with `RadialCharacteristics(m, ν, M)` and define the partition of actuators
1. actuator: b1
2. actuator: b1 
3. actuator: b2

The actuator setup on the left boundary is created with `setIOSetup!()`. 

```@example 2d_plate
config1  = RadialCharacteristics(1, 1, 4) # Set actuators characteristics b1, b2
config2  = RadialCharacteristics(0.9, 2, 20) # Set actuators characteristics b3

partition    = [1,2,3]
config_table = [config1, config1, config2]

setIOSetup!(actuation, plate, partition, config_table,  pos_actuators)
```

## Simulation
In this section, we define our heat conduction problem as an ordinary differential equation (ODE) that can be solved with `OrdinaryDiffEq.jl` or `DifferentialEquations.jl`. We specifiy the constant input signal $u_{1}(t)=u_{2}(t)=u_{3}(t)=2\cdot 10^{5}$ and develop function `heat_conduction!` as an interface for the ODE solver. 

```julia
function heat_conduction!(dθ, θ, param, t)
    u_in = 2e5 * ones(3)    # heat input
    diffusion!(dθ, θ, plate, property, boundary, actuation, u_in)
end
```

The inital temperature $\theta(0,x,y)=300$ Kelvin and the final simulation time is set to $300$ seconds. We use the numerical integration algorithm `KenCarp5` because it can handle stiff differential equations, and save the result for each second. 

```julia
using OrdinaryDiffEq
θinit = 300*ones(Nx*Ny) # Inital temperatures
tspan = (0.0, 300.0) # Time span
alg = KenCarp5()
prob = ODEProblem(heat_conduction!,θinit,tspan)
sol = solve(prob,alg, saveat=1.0)
```

## Evaluation
The actuator characteristics is displayed as a scatter plot.

```julia
xgrid = L/(2Nx) : L/Nx : L # Position in x-direction
ygrid = W/(2Ny) : W/Ny : W # Position in y-direction

using Plots
act_char, _ = getCharacteristics(actuation, pos_actuators)
scatter(xgrid, act_char, xlabel="Length L", ylabel="Spatial Characteristics",label=["Actuator 1" "Actuator 2" "Actuator 3"], legend=:bottomright, linewidth=3)
```

The evolution of the temperature is expressed by an animation.
```julia
anim = @animate for i in axes(sol,2)
    heatmap(xgrid, ygrid, reshape(sol[:,i], Nx, Ny)', xlabel="Length L", ylabel="Width W", title=string("Temperature at Time ",sol.t[i], " seconds") )
end
gif(anim, "plate_2d_animation.gif", fps=10)
```