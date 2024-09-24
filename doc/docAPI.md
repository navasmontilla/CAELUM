## API documentation

Here, we offer a summary of the main programming structures and functions for the user to implement additional functionalities. The implementation of other spatial reconstructions and Riemann solvers is straightforward.

### Data structures

In EHOW3D, the variables are organized using c structures. Below, there is a summary of the most relevant structures and variables (note that some variables are ommited in this document).

Cell data is contained in a structure of type **t_cell_**:

```c
struct t_cell_{
	int id;
	int l,m,n;//index in cartesian reference
	double *U; //this is the array of conserved variables. When using Euler: rho, rhou, rhov, E, rhophi
	double *U_aux; //this is an auxiliary array for RK stepping.
	double *Ue; //equilibrium state for atmospheric flow
	double *S; //source term
      double *S_corr; //correction of the source term
	double pres,prese,u_int; 
	double dx,dy,dz;
	double xc,yc,zc;	
	int n1,n2,n3,n4,n5,n6,n7,n8; //ID's of the nodes of the cell
	int w1_id,w2_id,w3_id,w4_id,w5_id,w6_id; //ID's of the walls of the cell
	...
	...	
      int st_sizeX, st_sizeY, st_sizeZ;	//stencil size
	int stX[9], stY[9], stZ[9];		//id's of the cells in the X and Y stencil.

};
```

Wall data is contained in a structure of type **t_wall_**:

```c
struct t_wall_{
	int id;
	int stencil; //this could be 1, 3, 5 or 7, depending on the method stencil
	double *UL, *UR; //array of reconstructed values on the left and right hand side of the wall, coming from (WENO/TENO) reconstruction
	double *fR_star,*fL_star; //array of numerical fluxes on the left and right hand side of the wall, provided by the Riemann solver
      double *ULe, *URe; //array of reconstructed values on the left and right hand side of the wall, coming from (WENO/TENO) reconstruction, for the EQUILIBRIUM
	double pRe,pLe; //equilibrium pressures
	int cellR_id, cellL_id; //id of the right and left cell
	t_cell *cellR, *cellL; //pointers to the left and right hand cells of the wall
	...
	
};
```

here, we would like to highlight the variables ```UL, UR``` as well as ```fR_star, fL_star```, which are the reconstructed conserved variables and the numerical fluxes on the left and right hand side of the wall.

Other global parameters are contained in the structures **t_mesh_** and **t_sim_**. There are additional structures to handle solid bodies.

### Numerical core 

The most relevant functions featuring the FV scheme are displayed below

```c
void update_cell(t_mesh *mesh, t_sim *sim);	//First order explicit Euler integration in time 
void update_cellK1(t_mesh *mesh, t_sim *sim);	//First step for SSPRK3
void update_cellK2(t_mesh *mesh, t_sim *sim);	//Second step for SSPRK3
void update_cellK3(t_mesh *mesh, t_sim *sim);	//Third step for SSPRK3
int compute_fluxes(t_mesh *mesh, t_sim *sim);	//Spatial reconstruction and computation of numerical fluxes
void update_solution(t_mesh *mesh, t_sim *sim, t_solid *solids, int rk_steps); //Core of the update in time
```

The main logic of the algorithm can be seen in **update_solution()**. It is displayed below (some parts have been omitted for the sake of clarity):

```c
void update_solution(t_mesh *mesh, t_sim *sim, t_solid *solids, int rk_steps){
	int k;
	for(k=1;k<=rk_steps;k++){
		if(k==1){
			compute_fluxes(mesh,sim);
			#if ST!=0&&EQUATION_SYSTEM==2
				compute_source(mesh);
			#endif
			update_dt(mesh,sim);
			if(rk_steps==1){
				update_cell(mesh,sim);
				(...)
			}else{
				update_cellK1(mesh,sim);
				(...)                               
			}
		}else if(k==2){
			compute_fluxes(mesh,sim);
			#if ST!=0&&EQUATION_SYSTEM==2
				compute_source(mesh);
			#endif
			update_cellK2(mesh,sim);
			(...)   				
		}else{ //k=3
			compute_fluxes(mesh,sim);
			#if ST!=0&&EQUATION_SYSTEM==2
				compute_source(mesh);
			#endif
			update_cellK3(mesh,sim);
			(...)     
		}
	}
}
```


### Spatial Reconstructions

The spatial reconstructions (on the left and right sides of an interface, L and R) are implemented in the functions below:

```c
double weno3L(double *phi);
double weno3R(double *phi);
double weno5L(double *phi);
double weno5R(double *phi);
double weno7L(double *phi);
double weno7R(double *phi); 
```

In spite of their name, each function features the WENO, TENO and optimal reconstructions. As an example, the 3-rd order reconstruction on the right hand side of a wall:

```c
double weno3R(double *phi){

	double g0, g1;		//gamma: optimal weight/	
	double w0, w1;		//omega: WENO weight
#if TYPE_REC < 2
	double b0, b1;		//beta: smoothness indicators
	double a0, a1;		//alpha
#if TYPE_REC == 1
	double c0, c1;
#endif
#endif
	double UR;

	g0=2.0/3.0;
	g1=1.0/3.0;

#if TYPE_REC == 0 //WENO

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));

	w0=a0/(a0+a1);
	w1=a1/(a0+a1);

#elif TYPE_REC == 1  //TENO

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=1.0/pow((b0+epsilon2),_Q_);
	a1=1.0/pow((b1+epsilon2),_Q_);

	c0 = a0/(a0 + a1);
      c1 = a1/(a0 + a1);

      c0 = c0 < _CT_ ? 0. : 1.;
      c1 = c1 < _CT_ ? 0. : 1.;

	a0 = g0*c0;
      a1 = g1*c1;

      w0 = a0/(a0 + a1);
      w1 = a1/(a0 + a1);

#else //Optimal

	w0=g0;
	w1=g1;
	
#endif

	UR=w0*(0.5*phi[1]+0.5*phi[0])+w1*(-0.5*phi[2]+1.5*phi[1]);
	return UR;
}
```

where ```phi[0],phi[1],phi[2]``` are the values of a conserved variable in each of the cells in the stencil.

### Riemann Solvers

The available Riemann solvers are given below.

For the **linear scalar equation**, we use an upwind flux definition, implemented in:
```c 
void compute_linear_flux(t_wall *wall,double *lambda_max);
```
For the **Burgers equation**, we also use an upwind flux definition, implemented in:
```c 
void compute_burgers_flux(t_wall *wall,double *lambda_max);
```
For **Euler equations**, the available solvers are:
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
```

The main logic of a solver function is to use the left and right reconstructed data (```wall->UL``` and ```wall->UR```) to compute the fluxes (```wall->fL_star``` and ```wall->fR_star```). These solvers are implemented using the x-split equations, that means that the vector variables must be rotated to the x-reference axis (```WL``` and ```WR```) before computing the numerical flux. Afterwards, the fluxes are rotated to the original direction (wall normal).


```c
/**Rotation of array of variables**/
WR[0]=wall->UR[0]; //Mass is not vectorial
// This is a simplification of the rotation matrix, only valid for cartesian mesh
WR[1]=wall->UR[1]*wall->nx+wall->UR[2]*wall->ny+wall->UR[3]*wall->nz;
WR[2]=-wall->UR[1]*wall->ny+wall->UR[2]*wall->nx+wall->UR[2]*wall->nz;
WR[3]=wall->UR[3]*wall->nx+wall->UR[3]*wall->ny-wall->UR[1]*wall->nz;
WR[4]=wall->UR[4]; //Energy is not vectorial
```

```c
/**Inverse rotation of the flux**/
wall->fR_star[0]=F_star[0]; //Mass is not vectorial
wall->fR_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny - F_star[3]*wall->nz;
wall->fR_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx + F_star[2]*wall->nz;
wall->fR_star[3]=F_star[3]*wall->nx + F_star[3]*wall->ny + F_star[1]*wall->nz;
wall->fR_star[4]=F_star[4]; //Energy is not vectorial
```