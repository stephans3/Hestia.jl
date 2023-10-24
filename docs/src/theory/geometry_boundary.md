# Geometry

Hestia is developed for simple geometries in one, two and three dimensions. At the moment these shapes are implemented:
- 1D rod: `HeatRod`
- 2D rectangle (or plate): `HeatPlate`
- 3D cuboid: `HeatCuboid`

A **HeatRod** has
- a length (in x-direction) and
- a number of discretization cells.

A **HeatPlate** has
- a length (in x-direction), its number of discretization cells (Nx) and
- a width (in y-direction), its disretization (Ny).

A **HeatCuboid** has
- a length (in x-direction), its number of discretization cells (Nx),
- a width (in y-direction), its disretization (Ny) and
- a height (in z-direction), its disretization (Nz).


## Boundary sides
The **one-dimensional rod** (HeatRod) has two boundary sides:
- `:west` and `:east`.

These boundary sides are single points.

The **two-dimensional plate** (HeatPlate) has four boundary sides:
- `:west`, `:east` and
- `:south`, `:north`.

These boundary sides are one-dimensional lines. 

The **three-dimensional cuboid** (HeatCuboid) has four boundary sides:
- `:west`, `:east`,
- `:south`, `:north` and
- `:underside`, `:topside`.

These boundary sides are two-dimensional areas. 
