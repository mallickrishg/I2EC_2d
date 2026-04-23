# I2EC: Integral Imposed Earthquake Cycles (2D)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs40623--026--02429--4-green)](https://doi.org/10.1186/s40623-026-02429-4)

**Authors:** Rishav Mallick & Sharadha Sathiakumar (Earth Observatory of Singapore)

**I2EC** is a MATLAB-based numerical modeling tool designed to simulate the spatiotemporal evolution of earthquake cycles in two-dimensional subduction zone geometries. It solves for mechanical equilibrium in the lithosphere-asthenosphere system to predict displacements, strains, and stresses. 

The codebase utilizes an extended Boundary Element Method (BEM) framework that reduces governing partial differential equations to a system of coupled ordinary differential equations (ODEs).

## Key Features

* **Mechanical Consistency Across Timescales:** Ensures consistency from individual earthquake cycles to geological periods by integrating long-term deformation rates with short-term cycle perturbations.
* **Support for Diverse Rheologies:** Accommodates spatially variable power-law viscous rheologies in the bulk mantle and rate-dependent friction on faults.
* **Advanced Numerical Schemes:**
    * **ODE Integration:** Solves the coupled system (IVP) using a 4th-order Runge-Kutta scheme for non-linear cases.
    * **Semi-Analytical Speedup:** Transforms linear viscoelastic problems into an eigenvalue problem, reducing computational time by orders of magnitude.
* **Optimized Stress Kernels:** Implements regularized stress kernels to eliminate unphysical instabilities common in standard hypersingular BEM formulations.
* **Geodynamical Integration:** Directly incorporates long-term strain rates derived from geodynamic solutions (e.g., Newtonian corner-flow approximations) as a loading term.

## Mathematical Overview

The code operates by splitting deformation into two components:
1.  **Long-term Component:** A quasi-static boundary value problem for background flow and steady fault slip.
2.  **Earthquake Cycle Component:** A time-dependent perturbation solved via an integral representation of the medium's elastic response.

The governing system of ODEs for the stress rate is:

$$\dot{\\sigma}_i(\mathbf{x}) = \int_{\Gamma} K_i(\mathbf{x};\xi)\cdot(v(\xi) - v^{\infty})d\xi + \iint_{\Omega} K_{ji}(\mathbf{x};\alpha)\cdot(\dot{\epsilon}_{\mathrm{v}_j}(\alpha) - \dot{\epsilon}_{\mathrm{v}_j}^{\infty}(\alpha))d\alpha$$

where $K$ represents the elastostatic kernels mapping fault slip and viscous strain to stress.

## Code Dependencies

This repository integrates and modifies components from these sources:
* **[MATLAB PDE Toolbox](https://www.mathworks.com/products/pde.html):** Used for 2D geometry discretization and mesh generation.
* **Barbot (2018):** Displacement and stress kernels for eigen sources in a linear elastic half-space (plane strain). [DOI: 10.1785/0120180058](https://doi.org/10.1785/0120180058)
* **[CutAndDisplace](https://github.com/Timmmdavis/CutAndDisplace):** Kernels for fault-related displacement and stress calculations.

## Citation

If you use this code in your research, please cite the original article:

> Mallick, R., & Sathiakumar, S. (2026). A fast numerical routine to model viscoelastic earthquake cycles in subduction zones. *Earth, Planets and Space*. [https://doi.org/10.1186/s40623-026-02429-4](https://doi.org/10.1186/s40623-026-02429-4)
