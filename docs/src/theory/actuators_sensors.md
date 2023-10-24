# [Actuators and Sensors](@id actuators)

Actuators are assumed only on the boundary sides of the geometrical object - not inside the object. It is possible to place more than one actuator on one boundary side. So, boundary sides with actuators are subdivided in **partitions** of equal size. In each partition we can specify wheter one actuator is assumed or not. It is not possible to specify more than one actuator per partition. If for instance 5 actuators are assumed on boundary `:south` then 5 partitions are created.

# Spatial Characterization

The n-th actuator has inside its partition the radially symmetric spatial characterization

$b_{n}(x) = m_{n} \exp(-M_{n} ~ \lVert x-x_{c}\rVert^{2~\nu_{n}})$

with
- scaling $m \in [0,1]$
- curvature matrix $M \in \mathbb{R}^3$,
- power $\nu \in \{1,2,3, \cdots\}$ and
- central point $x_{c}$ of the partition.

The central point of each partition is calculated internally, all other values have to be fixed. In many situations curvature matrix $M$ can be assumed as a scaled identity matrix, for example $M = 54 ~ I_{3 \times 3}$.

## Quick Overview
`HeatRod`s have single points as boundary sides: 
- only 1 actuator per side is possible
- no partitions
- spatial characterization $b$ is only a number 

`HeatPlate`s have four boundary sides as 1-dimensional lines:  
- more than 1 actuator per side is possible
- partitions are 1-dimensional intervals
- spatial characterization $b_{n}(x)$ is 1-dimensional function

`HeatCuboids`s have six boundary sides as 2-dimensional areas:  
- more than 1 actuator per side is possible
- partitions are 2-dimensional intervals
- spatial characterization $b_{n}(x)$ is 2-dimensional function, $x=(x_{1},x_{2})$


# Sensors

```@docs
measure
```