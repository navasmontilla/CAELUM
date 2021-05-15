# EHOW3D (Euler High Order WENO 3D)



## Introduction

3D compressible solver using 1-st, 3-rd, 5-th and 7-th order WENO, TENO and UWC reconstructions, SSPRK3 time updating and HLL/HLLC solvers. 

<figure style="text-align: center;">
  <img src="doc/panel.PNG" width="70%" alt="my alt text"/>
</figure>

This code allows the simulation of

- Linear scalar transport
- Burgers' equation
- Compressible Euler equations



## Quick start

### Compilation

To compile using *gcc* and the OpenMP *-fopenmp* flag, type:

```
$ gcc -lm -fopenmp euler3D.c -o euler3D
```

It can also be compiled using PGI NVIDIA compiler, with the OpenMP *-mp* flag:
```
$ pgcc -mp -o euler3D euler3D.c
```

### Configuration file

The configuration file required has the following format:

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


## Documentation

### Equations solved

There is the possibility of solving:

- Linear scalar transport, setting:
```c
#define LINEAR 1 
```
and defining the x, y and z velocities in the configuration file.
- Burgers equation, setting:
```c
#define BURGERS 1 
```
- Euler equations, setting:
```c
#define EULER 1 
```

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
### Solid domains

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

*This  feature is under construction yet.*


### Boundary conditions

The available boundary conditions are:

* 1: Periodic.

* 3: Transmissive. The numerical flux is set as the physical flux at the interface, using:
```c 
void compute_transmissive_euler(t_wall *wall, int wp)
```

* 4: Solid wall. Defined as a slip boundary condition which is based on the HLL flux, using:
```c 
void compute_solid_euler_hlle(t_wall *wall, double *lambda_max, int wp)
```




### Spatial reconstructions

Spatial reconstructions are implemented using 1D splitting. The available reconstructions are:

- UWC 3, 5 and 7
- WENO 3, 5 and 7
- TENO 5 (requires the selection of the CT constant!)

### Time integrator

The time stepping is done using a SSPRK3 method (3-rd order of accuracy).

### Riemann solvers

The available solvers are:

- HLL solver: 
```c 
void compute_euler_HLLE(t_wall *wall,double *lambda_max) 
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

### Parallel implementation

The code is parallelized using OpenMP. To set the number of threads, for example 24 threads, the following macro is used:
```c
#define NTHREADS 24
```
If compiling without the OMP flag, the code will run as a serial program.

## Authorship

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2018-2019 The authors.  

License type: Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Spain (CC BY-NC-ND 3.0 ES https://creativecommons.org/licenses/by-nc-nd/3.0/es/deed.en) under the following terms: 

- Attribution — You must give appropriate credit and provide a link to the license.
- NonCommercial — You may not use the material for commercial purposes.
- NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material unless explicit permission of the authors is provided. 

**Disclaimer:** This software is under development and it is distributed for research and/or academic purposes, WITHOUT ANY WARRANTY. In no event shall the authors be liable for any claim, damages or other liability, arising from, out of or in connection with the software or the use or other dealings in this software.

