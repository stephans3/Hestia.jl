# Getting Started with Hestia

In this tutorial you will learn how to simulate heat conduction phenomenas with **Hestia.jl**. Firstly, we take build a simulation for a 1-dimensional linear heat equation. And secondly, we show how to simulate 2-dimension *quasi-linear* heat conduction.

## Example: 1-dimensional linear heat conduction

In our first example, we show how to simulate the heat equation

```math
\frac{d}{dt} \theta(t,x) = \frac{\lambda}{c~\rho} \frac{d^2}{dx^2} \theta(t,x) \quad (t,x) \in (0,T) \times (0,L)
```

in a one-dimensional rod with length $L=0.2$ for $T=300$ seconds.  

#### Material
In this tutorial, we assume constant material properties
- thermal conductivity: $\lambda=50 \frac{W}{m ~ K}$,  
- specific heat capacity: $c=400 \frac{J}{kg ~ K}$ and
- mass density: $\rho=8000 \frac{kg}{m^3}$. 

We consider isotropic thermal conductivity because we model heat conduction only in one direction. 

```@example getting_started_1D
using Hestia
λ = 50    
ρ = 8000
c = 400 
prop = StaticIsotropic(λ, ρ, c)
```

!!! info
    Four types are offered to model the material properties: StaticIsotropic, DynamicIsotropic, StaticAnistropic, DynamicAnisotropic. More information about material properties can be found on page [Material and Physical Properties](theory/material_properties.md)


### Geometry

We consider a one-dimensional rod of length $L=0.2$ which is discretized by $40$ points. The one-dimensional rod is implemented in Hestia as `HeatRod` and it contains the length and number of grid points.

```@example getting_started_1D
L = 0.2    
Nx = 40  
rod = HeatRod(L, Nx)
```

#### Boundary sides

The one-dimensional rod has two boundary sides: on the left side and on the rigth side. The left side is denoted as `:west` and the right side as `:east`. On the boundary sides only Neumann-type boundary conditions can be used.
Here, we simulate our heat equation with thermal convection on boundary `:east` and thermal isolation on boundary side `:west`. The boundary condition with thermal convection is described mathematically as

```math
\lambda ~ \frac{d}{dx} \theta(t,x) = -h [\theta(t,x) - \theta_{amb}]
```

with `h` as heat transfer paramter and $\theta_{amb}$ as ambient temperature. We assume $h=10$ and $\theta_{amb}=300$.

```@example getting_started_1D
h = 10;
θamb = 300;
emission = Emission(h, 0, θamb)   # Only convection
boundary = Boundary(rod)
setEmission!(boundary, emission, :east); 
```

The thermal isolation on boundary side `:west` does not have to be specified in code explicitely. Internally, boundary side `:west` is initialized automatically as a zero-Neumann boundary condition.


!!! info
    Three geometric forms can be used: HeatRod for 1D, HeatPlate for 2D and HeatCuboid for 3D.  More information about the geometry and boundary sides can be found on page [Geometry and Boundaries](theory/geometry_boundary.md)

### Simulation
In the last step, we define our heat conduction problem as an ordinary differential equation (ODE) and solve it with `OrdinaryDiffEq.jl`. More information about how to solve differential equations can be found in the [DifferentialEquations.jl docs](https://docs.sciml.ai/DiffEqDocs/stable/). So we define the ODE function and call the Hestia function `diffusion!`. 

```@example getting_started_1D
function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, rod, prop, boundary)
end
```

We assume $\theta(0,x)=300 + 100 sin(\pi x/L)$ Kelvin as initial temperature and compute the simulation until $T=300$ seconds. Here, we use the `KenCarp5` numerical integrator because it is able to solve `stiff` differential equations. [Stiff equations](https://en.wikipedia.org/wiki/Stiff_equation) are explained on Wikipedia. 

```
using OrdinaryDiffEq
Ngrid = 1:1:Nx
θinit = 300 .+ 100sinpi.((Ngrid.-0.5)/Nx)
tspan = (0.0, 300.0)
alg = KenCarp5()
prob = ODEProblem(heat_conduction!,θinit,tspan)
sol = solve(prob,alg, saveat=1.0)
```

Finally, we plot our results with `Plots.jl` and compare the initial temperatures versus the final temperatures.  
```
using Plots
xgrid = L*(Ngrid.-0.5)/Nx
plot(xgrid, sol[:,1], title="Simulation of 1D heat conduction",xlabel="Length L in [m]", ylabel="Temperature in [K]", label="Initial")
plot!(xgrid, sol[:,end], label="Final")
```

## What is next?

Hestia is designed to simulate boundary controlled heat conduction phenomena in multiple dimensions (1D, 2D, 3D). On page [One-dimensional rod](tutorials/rod_1d.md) you will learn how to simulate 1D heat conduction including a heat supply on the left boundary side. On page [Two-dimensional plate](tutorials/plate_2d.md) we model quasi-linear heat conduction in a rectangle, in which the material properties depend on the temperature. 
