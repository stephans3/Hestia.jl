# Material

The material properties 
- volumetric mass **density** $\rho$,
- **specific heat capacity** $c$ and
- **thermal conductivity** $\lambda$

can be specified as 
- temperature-independent (here called *static*) or
- temperature-dependent (here called *dynamic*)

In case of **temperature-independent** material properties, the variables $\rho$, $c$ and $\lambda$ are defined as constant real values.

In case of **temperature-dependent** material properties, the variables $\rho$, $c$ and $\lambda$ are defined via `Vector`s. If the specific heat capacity is defined by 

$c(\theta) = 11 + 22~\theta + 33~\theta^2$

with temperature $\theta$, then the corresponding `Vector` is implemented as `c = [11, 22, 33]`. If at least one material property is temperature-dependent then the other properties have to be implemented as `Vector`s.

## Anisotropic heat conduction

Additional to the specification temperature-dependent vs. -independent the thermal conductivity can be assumed as **isotropic** or **anisotropic**. 

In case of *anisotropic* heat conduction the thermal conductivity of the geometrical object depends on the spatial direction. Mathematically noted, the thermal conductivity $\lambda$ is now a matrix or matrix-valued function, e.g. for cuboids

$\lambda = 
\begin{pmatrix}
\lambda_{x} & 0 & 0 \\
0 & \lambda_{y} & 0 \\
0 & 0 & \lambda_{z}
\end{pmatrix}$

or

$\lambda(\theta) =
\begin{pmatrix}
\lambda_{x}(\theta) & 0 & 0 \\
0 & \lambda_{y}(\theta) & 0 \\
0 & 0 & \lambda_{z}(\theta)
\end{pmatrix}.$

So for anisotropic heat conduction two (for 2D = plate) or three components (for 3D = cuboid) of thermal conductivity $\lambda$ have to be defined.

In conclusion, the material can be defined as
1. temperature-independent and isotropic: `StaticIsoProperty`
2. temperature-dependent and isotropic: `DynamicIsoProperty`
3. temperature-independent and anisotropic: `StaticAnisoProperty`
4. temperature-dependent and anisotropic: `DynamicAnisoProperty`