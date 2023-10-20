# Hestia.jl

<img src="https://raw.githubusercontent.com/stephans3/Hestia.jl/master/hestia_icon.svg" width="100">

Simulation of heat conduction problems in multiple dimensions with boundary control.



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7778449.svg)](https://doi.org/10.5281/zenodo.7778449)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stephans3.github.io/Hestia.jl/dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
<!--- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stephans3.github.io/Hestia.jl/stable) --->
<!---
[![Build Status](https://github.com/stephans3/Hestia.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stephans3/Hestia.jl/actions/workflows/CI.yml?query=branch%3Amain)--->
<!--- [![Coverage](https://codecov.io/gh/stephans3/Hestia.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/stephans3/Hestia.jl) --->

**Hestia.jl** is a [Julia](https://julialang.org) library to simulate controlled and uncontrolled heat conduction phenomena in multiple dimensions (1D, 2D, 3D). Hestia.jl offers many options to design your heat conduction simulation: geometry, material properties and specification of boundary sides.

**Geometries**:
* 1-dimensional rod of length 
* 2-dimensional plate of lenght L and width W
* 3-dimensional cuboid of length L, width W, height H


**Material Properties**:
* temperature-independent: **linear** heat conduction
* temperature-dependent: **quasi-linear** heat conduction
* isotropic
* anisotropic


**On boundary sides**:
* emit heat to environment via linear **heat transfer** and nonlinear **heat radiation**
* supply heat via actuators like heating elements
* measure surface temperature

The documentation will be updated soon.

For examples are listed in folder [**examples**](https://github.com/stephans3/Hestia.jl/tree/master/examples).

## Alternatives
If you look for general purpose approaches to simulate heat conduction, you may also take a look on
- [Trixi.jl](https://github.com/trixi-framework/Trixi.jl/)
- [Gridap.jl](https://github.com/gridap/Gridap.jl)
- [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl)


## Cooling-down processes: no actuation

- [**1D** with constant properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/cooling_down/1D/rod_iso_static_cool.jl): linear heat conduction
- [**1D** with temperature-dependent properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/cooling_down/1D/rod_iso_dynamic_cool.jl): quasi-linear heat conduction
- [**2D** with constant properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/cooling_down/2D/plate_iso_static_cool.jl) 
- [**2D** with temperature-dependent properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/cooling_down/2D/plate_iso_dynamic_cool.jl)
- [**3D** with constant properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/cooling_down/3D/cuboid_iso_static_cool.jl)
- [**3D** with temperature-dependent properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/cooling_down/3D/cuboid_iso_dynamic_cool.jl)


## Heating-up processes: with actuation

- [**1D** with constant properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/1D/rod_iso_static_heat.jl): linear heat conduction
- [**1D** with temperature-dependent properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/1D/rod_iso_dynamic_heat.jl): quasi-linear heat conduction
- [**2D** with constant properties (basic 1)](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/2D/plate_iso_static_heat_a1.jl) 
- [**2D** with constant properties (basic 2)](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/2D/plate_iso_static_heat_a2.jl) 
- [**2D** with constant properties (advanced 2)](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/2D/plate_iso_static_heat_b1.jl)
- [**2D** with constant properties (advanced 2)](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/2D/plate_iso_static_heat_b2.jl)
- [**2D** with temperature-dependent properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/2D/plate_iso_dynamic_heat.jl)
- [**3D** with constant properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/3D/cuboid_iso_static_heat.jl)
- [**3D** with temperature-dependent properties](https://github.com/stephans3/Hestia.jl/blob/master/examples/heating_up/3D/cuboid_iso_dynamic_heat.jl)


# How to cite
If you use Hestia.jl in your own research, please cite the following article:
```bibtex
@article{scholz2023hestia,
  title={Hestia. jl: A Julia Library for Heat Conduction Modeling with Boundary Actuation.},
  author={Scholz, Stephan and Berger, Lothar},
  journal={Simul. Notes Eur.},
  volume={33},
  number={1},
  pages={27--30},
  year={2023}
}
```

You can also cite our software on Zenodo:
```bibtex
@misc{scholz2023Hestia,
  title={Hestia.jl},
  author={Scholz, Stephan},
  year={2023},
  month={03},
  howpublished={\url{https://github.com/stephans3/Hestia.jl}},
  doi={10.5281/zenodo.7778449}
}
```