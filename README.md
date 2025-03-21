# MixedPrecisionIMC

Julia software package to solve thermal radiatiative transfer problems using the Implicit Monte Carlo method in a specified float-precision.

## Installation Instructions



## Dependencies
Random
Plots


This code uses an input deck system to allow users to configure many different options for simulations.

Example inputs are available in the src\inputs folder. These example files can be copied and have their keywords values modified to run problems.

A synopsis of available keywords to be used in input decks is provided here:

NAME: Enter a name for your simulation. Will be included on plot outputs.

PRECISION: Specify the float precision to use. Options are Float16, Float32, and Float64.

SEED: Specifify a seed value for the random number generator. Enter an integer.

GEOMETRY: 1D or 2D

MESHTYPE: UNIFORM or NONUNIFORM

If a uniform meshtype is selected specify
XSIZE: Width in x-dimension
YSIZE: Width of y-dimension (only used with 2D geometry)
DX: Cell width in x-dimension
DY: Cell width in y-dimension

If a nonuniform meshtype is selected specify
XMESHNODES: Array of nodal mesh points in x-dimension. Array with format \[x,x,x,x\] 
YMESHNODES: Array of nodal mesh points in y-dimension. Array with format \[x,x,x,x\]

Material properties
If using a 1D mesh, material regions are specified using an array with format [meshstart, mat1_end, ..., matN-1_end, meshend] 
If using a 2D mesh, material regions are specified using an array of tuples holding x and y-coordinates for rectangular regions [((mat1_xstart, mat1_xend),(mat1_ystart, mat1y_end)), ... ,((matN_xstart, matN_xend),(matN_ystart, matNy_end))]

SIGMA_A_REGS: Regions for each absorption cross-section
SIGMA_A_VALS: Constant values of absorption cross-section
SIGMA_A_POWERS: Exponent powers for temperature power law form of absorption cross-section. Set to zero for no temperature dependence.

SIGMA_S_REGS: Regions for each scattering cross-section
SIGMA_S_VALS: Constant values of scattering cross-section
SIGMA_S_POWERS: Exponent powers for temperature power law form of scattering cross-section. Set to zero for no temperature dependence.

BEE_REGS: Regions for heat capacity
BEE_VALS: Heat capacity values (Energy/Temperature)

RADSOURCE_REGS: Radiation source regions
RADSOURCE_VALS: Radiation source values

Boundary Conditions
Options are either REFLECT or VACUUM.
LEFTBC: Left Boundary Condition
RIGHTBC: Right Boundary Condition
TOPBC: Top Boundary Condition (Only used in 2D)
BOTTOMBC: Bottom Boundary Condition (Only used in 2D)

Timing
TIMESTEPPING: CONSTANT or RAMP
DT: Timestep size - Constant only
DT0: Initial timestep size - Ramp only
K: Timestep multiplication factor - Ramp only
DTMAX: Largest timestep size - Ramp only
ENDTIME: Simulation endtime

Particles
NINPUT: Number of particles to try to source in each step
NMAX: Maximum number of particles allowable in simulation

Temperature
Use same format as material region to find mesh body temperature distributions
T_INIT: Initial body temperature - A single value can be supplied to make temperature homogeneous throughout

Surface temperatures can be defined along segments of each boundary surface
In 1D: \[Tleft,Tright\], In 2D: [(BottomSeg1_End, ... ,BottomSegN_End), (TopSeg1_End, ... ,SegN_End), (LeftSeg1_End, ... ,LeftSegN_End), (RightSeg1_End, ..., RightSegN_End)]

T_SURFACE_REGS: Surface temperature segments 
T_SURFACE_VALS: Surface temperature values 

Constants - Change units to match problem
PHYS_C: Speed of light
PHYS_A: Radiation constant
ALPHA: Discretization parameter - typically leave as 1

Transport
RANDOMWALK: TRUE of FALSE - Use Randomwalk acceleration for optically thick media (1D Only currently)

Scaling
ENERGYSCALE: Scale factor for energy
DISTANCESCALE: Scale factor for distance

Outputs
PLOTVARS: RADENERGY, MATENERGY, TEMPERATURE, HEATMAP (2D Temperature) - Variable to plot
BENCHMARK: TRUE or FALSE
XBENCH: x-axis benchmark values with format: \[x,x,x,x\]
YBENCH: y-axis benchmark values with format: \[x,x,x,x\]
SAVEFIG: TRUE or FALSE - Save figure of final timestep value of plotted variable
SAVEANIMATION: TRUE or FALSE - Save animation of plotted variable


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://"simonbutson".github.io/MixedPrecisionIMC.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://"simonbutson".github.io/MixedPrecisionIMC.jl/dev/)
[![Build Status](https://github.com/"simonbutson"/MixedPrecisionIMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/"simonbutson"/MixedPrecisionIMC.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/"simonbutson"/MixedPrecisionIMC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/"simonbutson"/MixedPrecisionIMC.jl)

