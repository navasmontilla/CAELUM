# EHOW3D (Euler High Order WENO/TENO 3D)



## Introduction

This repository includes a very high-order CFD solver for academic purposes. The solver allows the use of 1-st, 3-rd, 5-th and 7-th order WENO, TENO and linear reconstruction in space. This code allows the simulation of:

- Linear transport
- Burgers' equation
- Compressible Euler equations 

Regarding the compressible Euler equations, the following features are allowed in EHOW3D:

- Compressible flows single and multicomponent (composed by a mixture of 2 gases)
- Gravitational source terms

The high order WENO and TENO schemes implemented in this code allow the simulation of turbulent flows using an Implicit Large Eddy Simulation (ILES) framework. ILES methods accurately reproduce the statistical behavior of turbulent flows. The  truncation errors of the scheme play the role of the common sub-grid scale filters used in traditional LES methods. High-fidelity simulations can be achieved when using this approach. A Kelvin-Helmholtz instability computed by EHOW-3D is shown below.

<figure style="text-align: center;">
  <img src="doc/panel.png" width="100%" alt="my alt text"/>
</figure>

## Installation

Clone the repository in your local computer:

```git clone https://github.com/navasmontilla/EHOW3D_public.git```

Compile the program as follows:

```make```

Note that the *Makefile* considers the following default definition 

```
DEBUG = 0
CFLAGS =  -Wall  -fopenmp
```

but the user can customize the compiling flags as desired. 

This software relies on other dependencies, listed below:

- [GCC](https://gcc.gnu.org/) or other C compiler
- [Python3](https://www.python.org/downloads/), for pre- and post-processing
- [Paraview](https://www.paraview.org/), for data visualization


## Quick guide

### Equations solved

There is the possibility of solving:

- Linear scalar transport, setting:
```c
#define LINEAR 1 
```
and defining the x, y and z velocities in the configuration file.

$$\frac{\partial u}{\partial t} + \lambda_x\frac{\partial u}{\partial x}+ \lambda_y\frac{\partial u}{\partial y}+ \lambda_z\frac{\partial u}{\partial z}= 0 $$

- Burgers equation, setting:
```c
#define BURGERS 1 
```

$$ \frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x}+ u\frac{\partial u}{\partial y}+ u\frac{\partial u}{\partial z}=0 $$

- Euler equations, setting:
```c
#define EULER 1 
```
<img src="https://latex.codecogs.com/svg.image?\dpi{150}&space;\frac{\partial}{\partial&space;t}\left[\begin{array}{c}&space;\rho\\&space;\rho&space;u\\&space;\rho&space;v\\&space;\rho&space;w\\&space;E&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;x}\left[\begin{array}{c}&space;\rho&space;u\\&space;\rho&space;u^{2}&plus;p\\&space;\rho&space;uv\\&space;\rho&space;uw\\&space;u(E+p)&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;y}\left[\begin{array}{c}&space;\rho&space;v\\&space;\rho&space;vu\\&space;\rho&space;v^{2}&plus;p\\&space;\rho&space;vw\\&space;v(E+p)&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;z}\left[\begin{array}{c}&space;\rho&space;w\\&space;\rho&space;wu\\&space;\rho&space;wv\\&space;\rho&space;w^{2}&plus;p\\&space;w(E+p)&space;\end{array}\right]=\left[\begin{array}{c}&space;0\\&space;0\\&space;-\rho&space;g\\&space;-\rho&space;g&space;w&space;\end{array}\right]" title="euler eqs" />

where

<img src="https://latex.codecogs.com/svg.image?\begin{matrix}&space;E=\rho\left&space;&space;(&space;\frac{1}{2}\mathbf{v}^2&space;&plus;&space;e&space;\right&space;),&space;\quad&space;p=(\gamma-1)\rho&space;e\equiv&space;(\gamma&space;-1)&space;\left&space;(&space;E-&space;\rho\frac{1}{2}\mathbf{v}^2&space;\right&space;)\\&space;\\\quad&space;H=\frac{E&plus;p}{\rho}\equiv&space;\frac{1}{2}\mathbf{v}^2&space;&plus;&space;h&space;,&space;\quad&space;h=e&plus;\frac{p}{\rho}\end{matrix}&space;" title="\begin{matrix} E=\rho\left ( \frac{1}{2}\mathbf{v}^2 + e \right ), \quad p=(\gamma-1)\rho e\equiv (\gamma -1) \left ( E- \rho\frac{1}{2}\mathbf{v}^2 \right )\\ \\\quad H=\frac{E+p}{\rho}\equiv \frac{1}{2}\mathbf{v}^2 + h , \quad h=e+\frac{p}{\rho}\end{matrix} " />


It is possible to run the two-component Euler equations, setting:
```c
#define MULTICOMPONENT 1
```
and set ```MULTI_TYPE=1``` to choose this Gamma formulation  $\phi =\gamma$ or  ```MULTI_TYPE=2``` to use this formulation $\phi =1/(\gamma-1)$  (see R. Abgrall, S. Karni, Computations of Compressible Multifluids, JCP 169 (2001)) for

$$\frac{\partial \rho\phi}{\partial t} + \frac{\partial \rho u\phi}{\partial x}=0$$


*It is possible to run a linear transport module within the euler equations, setting ```#define LINEAR_TRANSPORT 1 ```, but it is an old feature*




### The computational mesh

The computational mesh is constructed as follows:

<figure style="text-align: center;">
  <img src="doc/grid.png" width="60%" alt="my alt text"/>
</figure>

with cell numbers in green, wall numbers in red and node numbers in blue.

Each reference element (volume cell) is defined as follows:

<figure style="text-align: center;">
  <img src="doc/cell.png" width="30%" alt="my alt text"/>
</figure>

with wall numbers in red and node numbers in blue. These are defined in the structure ``` t_cell_``` as:
```c
struct t_cell_{
	//...
	int n1,n2,n3,n4,n5,n6,n7,n8; //ID's of the nodes 
	int w1_id,w2_id,w3_id,w4_id,w5_id,w6_id; //ID's of the walls 
	//...
}; 
```
### Solid domains (*under development*)

Solid domains are allowed by setting the macro: 
```c
#define ALLOW_SOLIDS 1 
```

and involve the following functions:

```c
int read_solids(t_mesh *mesh, t_solid *solids ); // Read STL triangulation files
int assign_cell_type(t_mesh *mesh,t_solid *solids); // Define ghost and solid cells
int update_stencils(t_mesh *mesh,t_sim *sim); // Update stencil sizes depending on ghost cell layers
int assign_wall_type(t_mesh *mesh); // Define calculation walls 
int assign_image_cells(t_mesh *mesh,t_solid *solids); // Define image points for immersed boundaries
int update_ghost_cells(t_sim *sim,t_mesh *mesh,t_solid *solids); // Update the values of ghost cells using image points
```

*This  feature is still under development.*


### Boundary conditions

The available boundary conditions are:

* 1: Periodic.

* 2: (not available)

* 3: Transmissive. The numerical flux is set as the physical flux at the interface, using:
```c 
void compute_transmissive_euler(t_wall *wall, int wp)
```

* 4: Solid wall. Defined as a slip boundary condition which is based on the HLL flux, using:
```c 
void compute_solid_euler_hlle(t_wall *wall, double *lambda_max, int wp)
```

### Initial conditions

Generally, the initial conditions can be implemented in:
```c 
int update_initial(t_mesh *mesh);
```

For Euler equations, the problem variables can be assigned for instance as follows:
```c
for(k=0;k<mesh->ncells;k++){
	p= ... ;
	rho= ... ;  
	phi= ... ;
	u= ... ;
	v= ... ;
	w= ... ;
	
	...
	
}
```

For Burgers' and linear transport, only one variable must be assigned:
```c
for(k=0;k<mesh->ncells;k++){
	cell[k].U[0]= ... ;
}
```


For Euler equations, **the initial conditions can also be set in the file** ```initial.out```:

```c 
VARIABLES = X, Y, Z, u, v, w, rho, p, phi 
CELLS = 40, 40, 40,
0.0075 0.0075 0.0075 0.0 0.0 0.0 1.0 1.0 0.0
0.0075 0.0075 0.0225 0.0 0.0 0.0 1.0 1.0 0.0
0.0075 0.0075 0.0375 0.0 0.0 0.0 1.0 1.0 0.0
...
```



### Spatial reconstructions

Spatial reconstructions are implemented using 1D splitting. The available reconstructions are:

- Linear 3, 5 and 7
- WENO 3, 5 and 7
- TENO 3, 5 and 7 (requires the selection of the CT constant!)

### Time integrator

The time stepping is done using a Strong Stability Preserving Runge-Kutta 3 (SSPRK3) method.

### Riemann solvers

The available solvers are:

- HLL solver: 
```c 
void compute_euler_HLLE(t_wall *wall,double *lambda_max) 
```
- HLLS solver: 
```c 
void compute_euler_HLLS(t_wall *wall,double *lambda_max) 
```
- HLLC solver: 
```c 
void compute_euler_HLLC(t_wall *wall,double *lambda_max)
```

**Note**: A positivity fix must be implemented in the HLLC solver to avoid stability issues

The x-split version of the solvers is implemented. It must be noted that the rotation matrices are simplified for the particular case of Cartesian geometries, leading to:
```c
	WL[1]=wall->UL[1]*wall->nx+wall->UL[2]*wall->ny+wall->UL[3]*wall->nz;
	WL[2]=-wall->UL[1]*wall->ny+wall->UL[2]*wall->nx+wall->UL[2]*wall->nz;
	WL[3]=wall->UL[3]*wall->nx+wall->UL[3]*wall->ny-wall->UL[1]*wall->nz;
```
where ``` WL``` is the 3D vector in the (cell) local coordinates and  ``` UL``` in the absolute coordinates.

### Input data

Input data must be placed in the folder ```case/```, and contains:

- *configure.input*: Configuration file.
- *initial.out* (optional): Initial condition file described above, same format and filetype than the ASCII *.out output file.
- *solid_list.txt* (optional): Contains the number of solid bodies and the path to their files.
- *solid1.txt* (optional): 1-st solid body file (ASCII STL format).
- ...
- *solidN.txt* (optional): N-th solid body file (ASCII STL format).


#### Configuration file

The configuration file *configure.input* has the following format:

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


### Output data

EHOW-3D allows printing data in *.vtk format and ASCII *.out files. To activate each of those output file types, use the macros:

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

The dumping time is set as *DumpTime* in the file *configure.input*. 

### Parallel implementation

The code is parallelized using OpenMP. To set the number of threads, for example 24 threads, the following macro is used:
```c
#define NTHREADS 24
```
If compiling without the OMP flag, the code will run as a serial program.


### Compilation

To compile using *gcc* and the OpenMP *-fopenmp* flag, type:

```
$ gcc -lm -fopenmp euler3D.c -o euler3D
```

It can also be compiled using PGI NVIDIA compiler, with the OpenMP *-mp* flag:
```
$ pgcc -mp -o euler3D euler3D.c
```



## Numerical results

Results are sorted in the following categories:

- [Benchmark #1: 1D Riemann problems](doc/benchmark1.md)
- [Benchmark #2: 2D Riemann problems](doc/benchmark2.md)
- [Benchmark #3: Two-fluid mixture problem](doc/benchmark3.md)
- [Benchmark #4: Convergence rate test](doc/benchmark4.md)
- [Benchmark #5: Taylor-Green vortex](doc/benchmark5.md)
- [Benchmark #6: Kelvin-Helmholtz instability](doc/benchmark6.md)
- [Benchmark #7: Colliding thermals](doc/benchmark7.md)
- Other results





## Authorship

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar 
 
Collaborators:
 - Javier Guallart Huertas (2022-present)

Copyright (C) 2019-2023 The authors and collaborators.  

License type: Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Spain (CC BY-NC-ND 3.0 ES https://creativecommons.org/licenses/by-nc-nd/3.0/es/deed.en) under the following terms: 

- Attribution — You must give appropriate credit and provide a link to the license.
- NonCommercial — You may not use the material for commercial purposes.
- NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material unless explicit permission of the authors is provided. 

If you want to contribute to this project or provide any feedback, please [contact us](mailto:anavas@unizar.es)! ;)

**Disclaimer:** This software is under development and it is distributed for research and/or academic purposes, WITHOUT ANY WARRANTY. In no event shall the authors be liable for any claim, damages or other liability, arising from, out of or in connection with the software or the use or other dealings in this software.

