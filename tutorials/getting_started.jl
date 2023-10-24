using Hestia 

#Material
λ = 50    # Thermal conductivity: constant
ρ = 8000  # Mass density: constant
c = 400   # Specific heat capacity: constant
prop = StaticIsotropic(λ, ρ, c)

# Geometry
L = 0.2     # Length
Nx = 40    # Number of elements: x direction
rod  = HeatRod(L, Nx)

### Boundaries ###
h = 10;
θamb = 300;
emission = Emission(h, 0, θamb)   # Only convection
boundary = Boundary(rod)
setEmission!(boundary, emission, :east)

### Simulation ###
function heat_conduction!(dθ, θ, param, t)
    diffusion!(dθ, θ, rod, prop, boundary)
end

Ngrid = 1:1:Nx
θinit = 300 .+ 100sinpi.((Ngrid.-0.5)/Nx)
tspan = (0.0, 300.0)

using OrdinaryDiffEq
alg = KenCarp5()
prob = ODEProblem(heat_conduction!,θinit,tspan)
sol = solve(prob,alg, saveat=1.0)

using Plots
xgrid = L*(Ngrid.-0.5)/Nx
plot(xgrid, sol[:,1], title="Simulation of 1D heat conduction",xlabel="Length L in [m]", ylabel="Temperature in [K]", label="Initial")
plot!(xgrid, sol[:,end], label="Final")
