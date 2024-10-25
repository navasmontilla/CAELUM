# CÆLUM: An academic high-order finite-volume solver for the compressible Euler equations and related models


## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation-instructions)
3. [Automated test](#automated-test)
4. [Example usage](#example-usage)
5. [Functionality documentation](#functionality-documentation)
6. [Other numerical results](#other-numerical-results)
7. [Community guidelines](#community-guidelines)
8. [Authorship](#authorship)


## Introduction

`CÆLUM` is a hydrodynamic solver for compressible flows, which solves the Euler equations of gas dynamics by means of very high order essentially non-oscillatory finite volume schemes in cartesian meshes. A key feature of this software is that it also allows the simulation of simplified dry atmospheric flows in the meso- and micro-scale. The core of the solver is written in C but pre-processing and post-processing can be easily done using Jupyter Notebook as provided in the examples. The Weighted Essentially Non-Oscillatory (WENO)  and Targeted Essentially Non-Oscillatory (TENO) spatial reconstructions, in combination with Runge-Kutta integrators, are used in `CÆLUM` to achieve a high resolution in smooth regions of the flow and to capture discontinuous solutions without spurious oscillations. These schemes also allow to adopt an Implicit Large Eddy Simulation (iLES) framework when computing turbulent flows. The numerical diffusion plays the role of the explicit dissipative sub-grid model, being this a favorable choice when seeking simplicity and computational efficiency. In summary, this hydrodynamic solver is applicable to turbulent compressible flows, as well as low Mach number flows such as simplified dry atmospheric flows and other wave propagation phenomena. 

It is designed from an academic perspective, where clarity and accessibility are prioritized. Therefore, it includes user-friendly pre-processing and post-processing tools based on Python and Jupyter Notebook. The repository comes with a series of Python scripts for the configuration and visualization of various example flows, ranging from simple scalar advection in 1D to more complex atmospheric or compressible flow cases in 3D. These scripts rely on standard libraries such as `numpy`, `matplotlib` and `pyvista`, the latter being a powerful module for data visualization and rendering. Additionally, several Jupyter Notebooks are included, where the steps for configuring the simulation tool and visualizing the results are explained in detail.

Regarding the implementation of the solver, `CÆLUM` uses OpenMP, a directive-based threading library that allows parallel computing in multi-core CPUs. The modification of the base code by the inlcusion of pragmas is minimal, preserving the readability and clarity of the code.

<figure style="text-align: center;">
  <img src="doc/panel.png" width="100%" alt="my alt text"/>
</figure>

## Installation instructions

Go to the desired location where you want to download CAELUM, e.g. ```me@myPc: SomeFolder/$``` and clone the repository in your local computer:

```git clone https://github.com/navasmontilla/CAELUM.git```

Go to CAELUM main folder ```me@myPc: SomeFolder/CAELUM/$``` and compile the program using *Makefile* as follows:

```make```

Compilation flags ```DEBUG``` and ```OMP```, for debugging and for multi-thread computing with OpenMP, are defined in ```Makefile```. ```DEBUG=0```  and ```OMP=1```  are set by default. For instance, if you seek serial computing, you can compile as follows:

```make OMP=0```

This software relies on other dependencies, listed below:

- [GCC](https://gcc.gnu.org/) or other C compiler
- [Python3](https://www.python.org/downloads/), for pre- and post-processing. The following packages need to be installed using ```pip install```:
	- *matplotlib*
	- *numpy*
	- *scipy*
 	- *imageio*
	- *pyvista*
	- *jupyter notebook*
	
To install Python3 and the above packages:

```
sudo apt update
sudo apt install -y python3 python3-pip
pip3 install pyvista matplotlib numpy scipy imageio notebook
```

Other possibility: you can also do it using a virtual environment in the main folder (```me@myPc: SomeFolder/CAELUM/$```):

```
sudo apt update
sudo apt install -y python3 python3-pip
python3 -m venv myenv
source myenv/bin/activate
pip3 install pyvista matplotlib numpy scipy imageio notebook
```

## Automated test

To check the functionality of the software, an automated test composed of 6 benchmarks can be run as follows from the main directory (```me@myPc: SomeFolder/CAELUM/$```):

```python3 python/autotest.py nt rec ord```

where 
- ```nt``` is the number of threads
- ```rec``` is the reconstruction method (0: WENO, 1: TENO, 2: Optimal)
- ```ord``` is the order of accuracy (**Only 1, 3, 5 and 7 are available**)

Example usage for 8 threads, using WENO and order 3: ```python3 python/autotest.py 8 0 3```
  
The benchmarks include:

- A convergence rate test for the linear scalar equation
- 4 Riemann Problems (RP) for the Euler equations
- The colliding thermals test case for the Euler equations with gravity

Within this test, the program is compiled and executed for every benchmark, giving a *Passed*/*Not Passed* output on the terminal after the execution. The results can be visualized in [autotest/autotest.md](autotest/autotest.md)

## Example usage

To get started in a user-friendly environment, some Jupyter Notebooks have been created:

-  [Getting started: setting up a linear transport case](python/caseScalar.ipynb)

-  [Computing and visualizing shock waves](python/caseRP3D.ipynb)

-  [First steps with atmospheric flows: the colliding thermals test case](python/caseCollidingBub.ipynb)

For a correct functionality, *Jupyter Notebook* must be launched from the software main directory, e.g. ```me@myPc: SomeFolder/CAELUM/$ jupyter notebook```

A more complete set of examples, scripted in Python, are can be found  [**here**](doc/docExampleCases.md). The Python scripts to generate the cases below must be launched from the software main directory (```me@myPc: SomeFolder/CAELUM/$```).


## Functionality documentation

For additional information about the equations solved, code organization and libraries, input and output files, etc., see the [the functionality documentation](doc/docFunctionality.md)

For a detailed documentation of the main programming structures and functions of the code, see [the API documentation](doc/docAPI.md)


## Other numerical results 

Some additional results are presented below:

- [Benchmark #1: Convergence rate test](doc/benchmark4.md)
- [Benchmark #2: Taylor-Green vortex](doc/benchmark5.md)
- [Benchmark #3: Kelvin-Helmholtz instability](doc/benchmark6.md)
- [Benchmark #4: Colliding thermals](doc/benchmark7.md)

## Community guidelines

### Opening an Issue

If you have encountered a bug or have a feature request, opening an issue is the first step.  Before opening a new issue, check the [Issues](https://github.com/navasmontilla/CAELUM/issues) section to see if your issue has already been reported. This helps prevent duplicates. You can create an issue by going [here](https://github.com/navasmontilla/CAELUM/issues/new/choose).

- For bugs: Include steps to reproduce the issue, the expected vs. actual behavior, error messages, etc.
- For feature requests: Describe the feature to be implemented.


### Submitting a Pull Request

Here is how to submit a pull request to suggest improvements or fix issues:

- Fork the base repository [here](https://github.com/navasmontilla/CAELUM/fork).
- Clone the forked repository, make changes, and push them back to the fork.
- Create a pull request between the base and forked repositories [here](https://github.com/navasmontilla/CAELUM/pulls).
- Wait for the pull request to be either approved or dismissed.


## Authorship

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

License type: The 3-Clause BSD License, under the following terms:

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

This software is provided by the copyright holders and contributors “as is” and any express or implied warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose are disclaimed. In no event shall the copyright holder or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.

If you want to contribute to this project or provide any feedback, please [contact us](mailto:anavas@unizar.es)! ;)

