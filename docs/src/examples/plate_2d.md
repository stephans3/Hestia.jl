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

$\lambda \frac{\partial \theta(t,x,y)}{\partial y} = -h (\theta(t,x,y) - \theta_{amb}) - \epsilon \sigma (\theta(t,x,y)^4 - \theta_{amb}^4)$

on boundary `:north` at $y=W$. The heat transfer coefficient $h$, emissivity $\epsilon$ and ambient temperature $\theta_{amb}$ are assumed to be equal for each boundary side. We assume an initial and ambient temperature of $300$ Kelvin. 

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

# 1. Material Properties and Geometry
In the first step we define the material properties and the geometry. As stated in the beginning we assume *isotropic*, temperature-*dependent* material properties and thus we create a `DynamicIsoProperty` with `createDynamicIsoProperty(λ, ρ, c)`. We create the geometry: here the two-dimensional plate `HeatPlate`.

```julia
λ = [20.0, 0.1]  # Thermal conductivity: temperature-dependend
ρ = [8000.0]     # Mass density: constant
c = [220.0, 0.6] # Specific heat capacity: temperature-dependend
property = createDynamicIsoProperty(λ, ρ, c)

L = 0.2    # Length
W = 0.1    # Width
Nx = 40    # Number of elements: x direction
Ny = 40    # Number of elements: y direction
plate    = HeatPlate(L, W, Nx, Ny, property)
```

# 2. Emission
Next, we define the emitted heat flux on the boundary sides `:west`, `:east` and `:north`. We assume the same linear heat transfer (convection) and nonlinear heat radiation for all boundaries as noted in the table above. 

```julia
h = 5.0       # Heat transfer coefficient
ϵ = 0.3       # Emissivity
θamb = 300.0; # Ambient temperature
emission  = createEmission(h, ϵ, θamb)   # Only convection
```
The emission has to be assigned for all boundary sides. We initialize the boundary with `initBoundary()` to yield a container which stores all emissions and set the emissions.

```julia
boundary_plate  = initBoundary(plate)
boundary_east = :east 
boundary_west = :west 
boundary_north = :north 
setEmission!(boundary_rod, emission, boundary_east)
setEmission!(boundary_rod, emission, boundary_west)
setEmission!(boundary_rod, emission, boundary_north)
```

# 3. Actuation
Firstly, we initialize `IOSetup` as a container for actuator and sensor setups. In this example, we assume three input signals on boundary side `:south`. 
```julia
plate_actuation = initIOSetup(plate)
pos_actuators = :south  # Position of actuators
```

The spatial characterization is defined by

$b(x) = m ~ \exp(-M ~ \lVert x - x_{c} \rVert^{2 \nu})$

with scaling $m \in [0,1]$, curvature matrix $M \in \mathbb{R}^{3 \times 3}$, power $\nu \in \{1,2,\cdots\}$ and central point $x_{c}$. Curvature matrix is a diagonal matrix and usually defined as a scaled identity. Central point $x_{c}$ is computed internally. The first both actuators (from left to right) have the same characterization $(m_{1}, M_{1}, \nu_{1})$ and the third actuator has characterization $(m_{2}, M_{2}, \nu_{2})$ as listed below


| Name               | Variable                 | Value / Formula                |
| :---               |    :----:                |  :---                          |
| Actuator b1, b2    | $(m_{1}, M_{1}, \nu_{1})$|  $(1.0,  4~I_{3\times3}, 3)$   |
| Actuator b3        | $(m_{2}, M_{2}, \nu_{2})$|  $(0.9, 20~I_{3\times3}, 2)$   |

```julia
# Actuator b1 and b2
m1 = 1.0;  # Scaling 
M1 = 4.0;  # Curvature 4*I 
ν1 = 3;    # Power

# Actuator b3
m2 = 0.9;   # Scaling 
M2 = 20.0;  # Curvature 20*I
ν2 = 2;     # Power
```

We define two characterizations with `setConfiguration(m, ν, M)` and define the partition of actuators
1. actuator: b1
2. actuator: b1 
3. actuator: b2

Finally, the actuator setup on the left boundary is created with `setIOSetup!()`. 

```julia
config1  = setConfiguration(m1, ν1, M1) # Set actuators characterization b1, b2
config2  = setConfiguration(m2, ν2, M2) # Set actuators characterization b3

partition    = [1,2,3]
config_table = [config1, config1, config2]

setIOSetup!(plate_actuation, plate, partition, config_table,  pos_actuators)
```
