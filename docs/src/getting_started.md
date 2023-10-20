# Getting Started with Hestia

In this tutorial you will learn how to simulate heat conduction phenomenas with **Hestia.jl**. Firstly, we take build a simulation for a 1-dimensional linear heat equation. And secondly, we show how to simulate 2-dimension *quasi-linear* heat conduction.

## 1. Example: 1-dimensional linear heat conduction

In our first example, we show how to simulate the heat equation

```math
\frac{d}{dt} \theta(t,x) = \frac{\lambda}{c~\rho} \frac{d^2}{dx^2} \theta(t,x) \quad (t,x) \in (0,T) \times (0,L)
```

in a one-dimensional rod with length $L=0.1$ for $T=100$ seconds.  First of all, we specify the material properties
- thermal conductivity: $\lambda$
- specific heat capacity: $c$
- mass density: $\rho$

as

```@example linHE1_properties
using Hestia
λ = 50    
ρ = 8000
c = 400 
property = createStaticIsoProperty(λ, ρ, c)
```

