## Functionality documentation

- [Code organization and libraries](#code-organization-and-libraries)
- [Configuration of the code for compilation](#configuration-of-the-code-for-compilation)
- [Equations solved](#equations-solved)
- [Input data](#input-data)
- [Output data](#output-data)
- [The computational mesh](#the-computational-mesh)
- [Boundary conditions](#boundary-conditions)
- [Spatial reconstructions](#spatial-reconstructions)
- [Time integrator](#time-integrator)
- [Riemann solvers](#riemann-solvers)

### Code organization and libraries

```
caelum/
├── main.c
├── Makefile
├── lib/
│   ├── closures.c
│   ├── closures.h
│   ├── definitions.h
│   ├── ibmutils.c
│   ├── ibmutils.h
│   ├── mathutils.c
│   ├── mathutils.h
│   ├── numcore.c
│   ├── numcore.h
│   ├── postproc.c
│   ├── postproc.h
│   ├── preproc.c
│   ├── preproc.h
│   ├── reconst.c
│   ├── reconst.h
│   ├── solvers.c
│   ├── solvers.h
│   └── structures.h
├── python/
│   ├── utils.py
│   ├── autotest.py
│   ├── caseExample.py
│   └── caseExample.ipynb
├── case/
│   ├── configure.input
│   ├── equilibrium.out
│   ├── initial.out
│   └── out/
└── README.md
```

#### main Directory
- **main.c**: Main source file for the CÆLUM solver.
- **Makefile**: Makefile for compiling the code.
- **lib/**: Directory containing library files and utilities.
- **python/**: Directory containing Python scripts for pre/post-processing and automation.
- **case/**: Directory for case configuration and output files.
- **README.md**: Documentation for the CÆLUM solver.

#### lib Directory
- **definitions.h**: Definitions and constants used across the code.
- **structures.h**: Data structures used in the code.
- **closures.c/h**: Functions and headers related to pressure closures.
- **ibmutils.c/h**: Utilities for immersed boundary method.
- **mathutils.c/h**: Mathematical utilities.
- **numcore.c/h**: Core numerical methods and routines.
- **postproc.c/h**: Post-processing utilities.
- **preproc.c/h**: Pre-processing utilities.
- **reconst.c/h**: Reconstruction methods for high-order schemes.
- **solvers.c/h**: Riemann solver routines.

#### python Directory
- **utils.py**: Utility functions.
- **autotest.py**: Script for automated testing of the solver.
- **caseExample.ipynb**: Script for generating case configurations.


### Configuration of the code for compilation

The file ```lib/definitions.h``` contains some definitions and constants that will be used for compilation. The most relevant for the user are listed below. Here, we provide a summary with links to other sections where some of these features are explained more in detail.

#### Reconstruction Method

```c
#define TYPE_REC 0
```

- *Description*: Defines the type of high order reconstruction method used. See more info [here](#spatial-reconstructions).
- *Possible Values*:
  - `0`: WENO (Weighted Essentially Non-Oscillatory).
  - `1`: TENO (Targeted Essentially Non-Oscillatory).
  - `2`: Optimal Polynomial Reconstruction.

#### Equation System

```c
#define EQUATION_SYSTEM 2
```

- *Description*: Specifies the type of equation system to solve. See more info [here](#equations-solved).
- *Possible Values*:
  - `0`: Linear advection equation.
  - `1`: Burgers' equation.
  - `2`: Compressible Euler equations.

#### Other features for Euler equations

```c
#define ST 3
```

- *Description*: Controls the inclusion of source terms for the Euler equations.  See more info [here](#equations-solved).
- *Possible Values*:
  - `0`: Source terms OFF.
  - `1`: Source terms ON (augmented version, needs the use of HLLS).
  - `2`: Source terms ON (perturbation version, needs the use of HLL).
  - `3`: Source terms ON (perturbation version and total energy is conserved, needs the use of HLL).


```c
#define MULTICOMPONENT 0
```

- *Description*: Indicates whether single or multicomponent Euler equations are used. See more info [here](#equations-solved).
- *Possible Values*:
  - `0`: Single component Euler equations.
  - `1`: Multicomponent Euler equations (two components with different specific heat ratio $\gamma$).


```c
#define MULTI_TYPE 2
```

- *Description*: Defines the type of multicomponent formulation used. See more info [here](#equations-solved).
- *Possible Values*:
  - `1`: Gamma ($\gamma$) formulation.
  - `2`: $\frac{1}{\gamma - 1}$ formulation (recommended, as per R. Abgrall and S. Karni in JCP 169 (2001)).
  
  
```c
#define ALLOW_SOLIDS 0
```

- *Description*: Defines whether or not solid objects are allowed).
- *Possible Values*:
  - `0`: No solid cells allowed.
  - `1`: Solid objects are defined via *stl* files (not recommended; only for very simple geometries).
  - `2`: Solid cells can be defined through a list as an input file.
  - `3`: Solid objects are defined using Level-Set approach, using the signed distance function (SDF) as input.

#### Solver Selection for Euler equations

```c
#define SOLVER 0
```

- *Description*: Specifies the solver used for the numerical method. See more info [here](#riemann-solvers).
- *Possible Values*:
  - `0`: HLL (Harten-Lax-van Leer) solver.
  - `1`: HLLC (Harten-Lax-van Leer Contact) solver.
  - `2`: HLLS (Harten-Lax-van Leer Source) solver.

#### OpenMP Configuration

```c
#define NTHREADS 4
```

- *Description*: Defines the number of OpenMP threads used for parallel execution.
- *Possible Values*: Any integer value representing the number of threads. For example:
  - `4`: Use 4 threads for parallel computation.

#### Output Files Configuration


```c
#define WRITE_VTK 1
```

- *Description*: Controls whether a VTK file is generated for visualization.  See more info [here](#output-data).
- *Possible Values*:
  - `0`: VTK file is not generated.
  - `1`: VTK file is generated.


```c
#define WRITE_LIST 1
```

- *Description*: Controls whether a list (`*.out`) file is generated. See more info [here](#output-data).
- *Possible Values*:
  - `0`: List file is not generated.
  - `1`: List file is generated.

#### Variable Printing (VTK)

- *Description*: Defines which variables are printed to VTK files for visualization.

```c
#define print_RHO 0
#define print_VELOCITY 1
#define print_ENERGY 0
#define print_PRESSURE 0
#define print_OVERPRESSURE 1
#define print_SOLUTES 0
#define print_POTENTIALTEM 1
```

- *Possible Values for Each Variable*:
  - `0`: Do not print this variable.
  - `1`: Print this variable.

- *Variables*:
  - `print_RHO`: Density.
  - `print_VELOCITY`: Velocity vector.
  - `print_ENERGY`: Energy.
  - `print_PRESSURE`: Pressure.
  - `print_OVERPRESSURE`: Pressure difference from the equilibrium pressure.
  - `print_SOLUTES`: Passive solute concentration.
  - `print_POTENTIALTEM`: Potential temperature.

#### Initial Data Input

```c
#define READ_INITIAL 1
```

- *Description*: Controls how the initial data is loaded. See more info [here](#input-data).
- *Possible Values*:
  - `1`: Initial data is read from a file.
  - `2`: Initial data is set programmatically in the function `update_initial()`.
  
#### Physical constants

```c
#define PI 3.141592653589793
#define _g_ 9.8
#define _gamma_ 1.4
#define _R_ 287.058
#define _p0_ 1.0E5
```

- *Variables*:
  - `PI`: Pi number.
  - `_g_`: Gravity acceleration.
  - `_gamma_`: Specific heats ratio.
  - `_R_`: Ideal gas constant.
  - `_p0_`: Reference atmospheric pressure.


### Equations solved

There is the possibility of solving:

#### Linear scalar transport

To run this model, we must define:

```c
#define EQUATION_SYSTEM 0 
```

and the following equation is solved:

$$\frac{\partial u}{\partial t} + v_x\frac{\partial u}{\partial x}+ v_y\frac{\partial u}{\partial y}+ v_z\frac{\partial u}{\partial z}= 0 $$

where $v_x$, $v_y$ and $v_z$ the are x, y and z velocities.

#### Burgers equation

To run this model, we must define:

```c
#define EQUATION_SYSTEM 1 
```

and the following equation is solved:

$$ \frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x}+ u\frac{\partial u}{\partial y}+ u\frac{\partial u}{\partial z}=0 $$

#### Euler equations:
To run this model, we must define:

```c
#define EQUATION_SYSTEM 2 
```

and the following system of equations is solved:

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0 \quad\hbox{Continuity}$$

$$\frac{\partial (\rho \mathbf{v})}{\partial t} + \nabla \cdot \left(\rho \mathbf{v} \otimes \mathbf{v} + p \mathbf{I}\right) = \rho \mathbf{g} \quad\hbox{Momentum}$$

$$\frac{\partial E}{\partial t} + \nabla \cdot \left((E + p) \mathbf{v}\right) = \rho \mathbf{v} \cdot \mathbf{g} \quad\hbox{Energy}$$



where $\rho$ is density, $\mathbf{v}$ is the velocity vector, $p$ is pressure and $\mathbf{g}=(0,0,g)^T$ is the gravitational acceleration vector. The energy is defined as  the sum of kinetic and internal energy

$$E=\rho(\frac{1}{2}\mathbf{v}^2+e)$$

One should note the relations $p=(\gamma-1)\rho e\equiv (\gamma-1)(E-\frac{1}{2}\rho\mathbf{v})$.

- When setting ```#define ST 0```, we assume $g=0$. 

- When setting  ```#define ST 1```, we consider non-zero gravity and need to use the solver HLLS. 

- When setting  ```#define ST 2```, we consider non-zero gravity and need to use the solver HLL in fluctuation version. 

- When setting ```#define ST 3```  we consider non-zero gravity and need to use the solver HLL in fluctuation version. Besides,the equation for the conservation of energy is solved in fully conservative form, defining energy as $E_T=\rho(\frac{1}{2}\mathbf{v}+e+gz)$, yielding to

$$\begin{align}
\frac{\partial E_T}{\partial t} + \nabla \cdot \left((E_T + p) \mathbf{v}\right) &= 0 
\end{align}$$


It is possible to run the two-component Euler equations, setting:
```c
#define MULTICOMPONENT 1
```
and set ```MULTI_TYPE=1``` to choose this Gamma formulation  $\phi =\gamma$ or  ```MULTI_TYPE=2``` to use this formulation $\phi =1/(\gamma-1)$  (see R. Abgrall, S. Karni, Computations of Compressible Multifluids, JCP 169 (2001)) for

$$\frac{\partial \rho\phi}{\partial t} + \frac{\partial \rho u\phi}{\partial x}=0$$


### Input data

Inside the **case/** directory we will find all files corresponding to the simulation case, we need as input files:

- **configure.input**: Input configuration file for running a case, that includes the global configuration and is of the following form:

```
/////SIMULATION_SETUP//////
FinalTime		0.2
DumpTime		0.05
CFL			0.25
Order			5

////////MESH_SETUP/////////
xcells			80
ycells			100
zcells			80
SizeX			0.80
SizeY			1.0
SizeZ			0.80

///////BOUNDARY_COND///////
Face_1(-y)			3
Face_2(+x)			3
Face_3(+y)			3
Face_4(-x)			3
Face_5(-z)			3
Face_6(+z)			3

///////LINEAR_TRANSPORT///////(if_applicable)
u_x                     1.0
u_y                     1.0
u_z                     1.0
```

- **initial.out**: Input file for initial conditions, that is of the following form (for scalar and Euler equations respectively):
``` 
VARIABLES = X, Y, Z, u 
CELLS = 40, 40, 40,
0.0075 0.0075 0.0075 0.0 
0.0075 0.0075 0.0225 0.0
0.0075 0.0075 0.0375 0.0 
...
```
 
``` 
VARIABLES = X, Y, Z, u, v, w, rho, p, phi 
CELLS = 40, 40, 40,
0.0075 0.0075 0.0075 0.0 0.0 0.0 1.0 1.0 0.0
0.0075 0.0075 0.0225 0.0 0.0 0.0 1.0 1.0 0.0
0.0075 0.0075 0.0375 0.0 0.0 0.0 1.0 1.0 0.0
...
 ```
The initial conditions can also be set in, when setting ```READ_INITIAL = 0```:
```c 
int update_initial(t_mesh *mesh);
```

- **equilibrium.out**: Input file for equilibrium state (only when considering atmospheric cases). Similar structure than above.


- **solid_cells.input**: Input file for the definition of solid cells (`0` if solid) if `#define ALLOW_SOLIDS 2`
``` 
VARIABLES = X, Y, Z, sld 
CELLS = 40, 40, 40,
0.0075 0.0075 0.0075 0 
0.0075 0.0075 0.0225 1
0.0075 0.0075 0.0375 0
...
```

### Output data

This software allows printing data in *.vtk format and ASCII *.out files. To activate each of those output file types, use the macros:

```c
#define WRITE_VTK 1  //print *.vtk
#define WRITE_LIST 1 //print ASCII *.out
```

For *.vkt files, it is posible to choose the variables to print by means of additional macros. For instance, if we want to print X,Y,Z momentum and pressure, do:
```c
#define print_RHO 0
#define print_MOMENTUM 1
#define print_ENERGY 0
#define print_PRESSURE 1
#define print_OVERPRESSURE 0
#define print_SOLUTES 0
#define print_POTENTIALTEM 0
```

The time lapse for writing files is set as *DumpTime* in the file *configure.input*. 

### The computational mesh

The computational mesh is Cartesian and is constructed as follows:

<figure style="text-align: center;">
  <img src="../doc/grid.png" width="60%" alt="my alt text"/>
</figure>

with cell numbers in green, wall numbers in red and node numbers in blue.

Each reference element (volume cell) is defined as follows:

<figure style="text-align: center;">
  <img src="../doc/cell.png" width="30%" alt="my alt text"/>
</figure>

with wall numbers in red and node numbers in blue. These are defined in the corresponding data structures in *structures.h* (see [the API documentation](docAPI.md)).

### Boundary conditions

The available boundary conditions are:

* 1: Periodic
* 2: User defined
* 3: Transmissive (Euler). The numerical flux is set as the physical flux at the interface.
* 4: Solid wall (Euler). Defined as a slip boundary condition which is based on the HLL flux.

They are set by using the numbers above in the configuration file **configure.input**.

### Spatial reconstructions

Spatial reconstructions are implemented using 1D splitting. The available reconstructions are:

- Linear 3, 5 and 7
- WENO 3, 5 and 7
- TENO 3, 5 and 7

To select the spatial reconstruction method, use:
 ```c
#define TYPE_REC 0 //This is 0 for WENO, 1 for TENO and 2 for UWC
```

Note that **only orders 1, 3, 5 and 7 are available**.

### Time integrator

The time stepping is done using a Strong Stability Preserving Runge-Kutta 3 (SSPRK3) method when the spatial order is greater than 1, or with a 1-st order explicit Euler method, when the spatial order is 1.

### Riemann solvers

For the **linear scalar equation**, we use an upwind flux definition.

For the **Burgers equation**, we also use an upwind flux definition.

For **Euler equations**, the available solvers are:
- HLL solver
- HLLS solver
- HLLC solver

and can be configured as follows:
```c
#define SOLVER 0 //0: HLL solver, 1: HLLC solver, 2: HLLS solver
```
