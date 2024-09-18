---
title: 'EHOW3D: An academic High Order Finite Volume solver for hyperbolic problems'
tags:
  - fluid mechanics
  - compressible flow
  - atmospheric flow
  - finite volume
  - high order
authors:
  - name: Adrian Navas-Montilla
    orcid: 0000-0002-3465-6898
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Isabel Echeverribar
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Fluid Dynamics Technologies, i3A-University of Zaragoza, Spain
   index: 1
 - name: Institution Name, Country
   index: 2
date: 4 August 2024
bibliography: paper.bib

---

# Summary

In the last decade, very high order numerical methods have become very popular within the Computational Fluid Dynamics (CFD) and Numerical Weather Prediction (NWP) communities. Traditional low order --i.e. equal or lower than 2-- methods have been widely adopted by industry and academia to simulate fluid flows, but they show important limitations as they feature high diffusion and dispersion errors which may lead to non-physical solutions [@FERRER2023108700;@wang2013high]. On the other hand, high order schemes are able to provide accurate solutions, with small diffusion and dispersion errors, and feature a high computational efficiency. They are well suited for multi-scale problems, e.g. turbulent flows, flows involving sharp gradients such as shocks, acoustic phenomena, etc. High order solvers  have gained attention due to their notable performance when running on modern High-Performance Computing (HPC) architectures.


`EHOW3D` is a hydrodynamic solver for compressible flows, which solves the Euler equations of gas dynamics by means of very high order essentially non-oscillatory finite volume schemes in cartesian meshes. The Weighted Essentially Non-Oscillatory (WENO) and Targeted Essentially Non-Oscillatory (TENO) spatial reconstructions, in combination with Runge-Kutta integrators, are used in `EHOW3D` to achieve a high resolution in smooth regions of the flow and to capture discontinuous solutions without spurious oscillations. These schemes also allow to adopt an Implicit Large Eddy Simulation (iLES) framework when computing turbulent flows [@grinstein2007implicit]. The numerical diffusion plays the role of the explicit dissipative sub-grid model, being this a favorable choice when seeking simplicity and computational efficiency. In summary, this hydrodynamic solver is applicable to compressible turbulence, shock waves, and shock-induced turbulence, as well as low Mach number flows such as simplified dry atmospheric flows and other wave propagation phenomena.

 

# Statement of need

Given the limitations of traditional low-order numerical methods in accurately simulating complex fluid flows —especially in scenarios involving multi-scale phenomena like turbulence and shocks— there is a critical need for advanced, high-order computational tools. The increasing complexity of engineering and scientific problems, coupled with the availability of modern HPC architectures, necessitates a shift towards methods that can deliver both high accuracy and computational efficiency. It is thus required to popularize such methods and make them available for the scientific community, in particular, for students and novel researchers who are often limited by the accessibility of advanced tools.

`EHOW3D` is designed to be used both by researchers and by students in the field of CFD, atmospheric flows and related disciplines. In fact, the simplicity and compactness of the code make it particularly suitable for those taking their first steps in the numerical simulation of flows using high-order methods. It is designed from an academic perspective, where clarity and accessibility are prioritized. Therefore, it includes user-friendly pre-processing and post-processing tools based on Python and Jupyter Notebook. The repository comes with a series of Python scripts for the configuration and visualization of various example flows, ranging from simple scalar advection in 1D to more complex atmospheric or compressible flow cases in 3D. These scripts rely on standard libraries such as `numpy`, `matplotlib` and `pyvista`, the latter being a powerful module for data visualization and rendering. Additionally, several Jupyter Notebooks are included, where the steps for configuring the simulation tool and visualizing the results are explained in detail.

A key feature of this software is that it also allows the simulation of simplified dry atmospheric flows in the meso- and micro-scale by means of the compressible Euler equations with source term. We can find in the literature various formulations of the Euler equations for the construction of operational and research NWP dynamical cores. The most common version of the Euler equations for the NWP community is given by the equation for the conservation of mass, momentum and potential temperature [@giraldo2008study]. The latter equation is derived from the equation for the conservation of energy, assuming the flow is adiabatic. In contrast, `EHOW3D` uses the original Euler equations, composed by the equations for the conservation of mass, momentum and energy. This approach offers several benefits, such as the ability to conserve mass and energy with machine accuracy when using a suitable discretization. Additionally, many numerical advances developed by the CFD community, such as spatial reconstruction schemes and Riemann solvers, can be easily adapted for NWP [@giraldo2008study].

Regarding numerics, the simulation code uses WENO [@jiang1996efficient], TENO [@fu2019very], and optimal polynomial reconstructions on Cartesian meshes. For simplicity, we follow the strategy outlined in @zhang2011, which is based on the midpoint rule and employs independent 1D reconstructions in each Cartesian direction. This approach allows us to avoid performing multi-dimensional reconstructions and Gaussian integration at the cell faces, thereby drastically reducing the complexity of the algorithms as well as computational expenses. This is done at the cost of not achieving a genuinely high order of accuracy, which is not critical when computing shocked problems, underresolved flows, or flows with discontinuities and sharp gradients [@zhang2011]. As a result, we provide a simple and versatile computational code that can be applied to a wide variety of problems and enables iLES.

# The model equations

This simulation code considers the compressible Euler equations with gravitational source term, given by
\begin{align}
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) &= 0 \tag{Continuity} \\
\frac{\partial (\rho \mathbf{v})}{\partial t} + \nabla \cdot \left(\rho \mathbf{v} \otimes \mathbf{v} + p \mathbf{I}\right) &= \rho \mathbf{g} \tag{Momentum} \\
\frac{\partial E}{\partial t} + \nabla \cdot \left((E + p) \mathbf{v}\right) &= \rho \mathbf{v} \cdot \mathbf{g} \tag{Energy}
\end{align}
where $\rho$ is density, $\mathbf{v}$ is the velocity vector, $p$ is pressure and $\mathbf{g}=(0,0,g)^T$ is the specific gravity force. The energy is defined as  the sum of kinetic and internal energy $E=\rho(\frac{1}{2}\mathbf{v}+e)$. 

Additionally, scalar transport can also be considered as the following equation is implemented
$$\frac{\partial u}{\partial t} + \nabla \cdot ( \mathbf{v} u) = 0$$
where $u$ is the transported quantity and $\mathbf{v}$ is the advection velocity.  It is also possible to compute the Burgers equation, i.e. $\mathbf{v}=1/2(u,u,u)^T$. 

Further details on the model equations and numerical resolution methods can be found in [@NAVASMONTILLA2023; @NAVASMONTILLA2024].

# Previous work using the software

`EHOW3D` was developed and first used in @NAVASMONTILLA2023 where we presented a novel methodology to construct very high order well-balanced schemes for the computation of the Euler equations with gravitational source term, with application to numerical weather prediction (NWP). The proposed approach uses augmented Riemann solvers, which allow for the preservation of the exact equilibrium between fluxes and source terms at cell interfaces. In particular, the augmented HLL solver (HLLS) was considered. The objective of this work was twofold: first, to assess the use of augmented solvers in computing the Euler equations with gravity, and second, to evaluate the performance of the novel TENO reconstruction for NWP.

Afterwards, a more traditional implementation of the Euler equations with gravity based on the use of fluctuation variables --i.e. avoiding particular well-balancing techniques such as augmented solvers-- was implemented in different forms to also allow for the conservation of total energy. The resulting schemes were evaluated in terms of their spectral resolution for the computation of turbulent flows in the atmosphere. The results were presented in @NAVASMONTILLA2024, with the objective of determining whether or not these models can be used to build an implicit Large Eddy Simulation framework, shedding light on their potential advantages or limitations in representing under-resolved atmospheric processes in the meso- and micro-scales. 


# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

Part of this work has been funded by the Spanish Ministry of Science, Innovation and Universities - Agencia Estatal de Investigación (10.13039/501100011033) and FEDER-EU under project-nr. PID2022-141051NA-I00.

# References