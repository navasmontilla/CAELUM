/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - structures.h

Content:
  -This header file contains the structures used in the code.

*/

#ifndef STRUCTURES_H
#define STRUCTURES_H
#include "definitions.h"


////////////////////////////////////////////////////
//////////////  S T R U C T U R E S  ///////////////
////////////////////////////////////////////////////


typedef struct t_node_ t_node;
typedef struct t_cell_ t_cell;
typedef struct t_wall_ t_wall;
typedef struct t_mesh_ t_mesh;
typedef struct t_sim_ t_sim;
typedef struct t_solid_ t_solid;
typedef struct t_stl_ t_stl;
typedef struct t_triangle_ t_triangle;


struct t_node_{
	int id;
	double x,y,z; //coordinates of the node

};


struct t_cell_{
	int id;
	int l,m,n;//index in cartesian reference
	double *U; //array of conserved variables. When using Euler: rho, rhou, rhov, E, rhophi
	double *U_aux; //aux array of conserved variables for RK stepping
	double *Ue; //array of conserved variables for the equilibrium state
	double *S; //array of source terms
	double S_corr; //source term correction for well-balancing
	double pres,prese,u_int,pmax; 
	double dx,dy,dz; //cell sizes
	double xc,yc,zc; //cell centers
	
	int n1,n2,n3,n4,n5,n6,n7,n8; //ID's of the nodes of the cell
	int w1_id,w2_id,w3_id,w4_id,w5_id,w6_id; //ID's of the walls of the cell
	t_wall *w1, *w2, *w3, *w4, *w5, *w6; //pointer to cell walls

	int type; //Cell type: 0=solid, 1=fluid, 2=ghost
      double sdf,nr[3]; //signed distance to solids
	double xim,yim,zim;
	int ni[8]; //id: for each ghost cell, neighbors of image points.
	double li[8],distabs;
	int distsolx,distsoly,distsolz; //distance to solid cells in cell units
	int solid_id,triangle_id;
	int out; //auxiliary variable for assing_cell_type() . It indicates that a cell is outside of a surface
	t_triangle *tri;
	int st_sizeX, st_sizeY, st_sizeZ;	//stencil size
	int stX[9], stY[9], stZ[9];		//id's of the cells in the X and Y stencil.
      int troubled;

};


struct t_wall_{
	int id;
	int stencil; //this could be 1, 3, 5 or 7, depending on the method stencil
	double *UL, *UR; //array of reconstructed values on the left and right hand side of the wall, coming from (WENO/TENO) reconstruction
	double *fR_star,*fL_star; //array of numerical fluxes
	double *ULe, *URe; //array of reconstructed equilibrium values on the left and right hand side of the wall, coming from (WENO/TENO) reconstruction, for the EQUILIBRIUM
	double pRe,pLe; //equilibrium pressures on the left and right hand side of the wall
	int cellR_id, cellL_id; //id of the right and left cell
	t_cell *cellR, *cellL; //pointers to the left and right hand cells of the wall
	double nx, ny, nz; //wall normals
	double z; //height in z direction
      int wtype, boundId; //wtype: 1 for inner walls, 3 for transmissive boundary walls and 4 for solid walls
      double vel; //advection velocity projected onto the wall normal
      int boundary; //=1 if wall belongs to a boundary
};



struct t_mesh_{
	int xcells, ycells, zcells; 
	double dx, dy, dz;
	double Lx, Ly, Lz;
	double u_x, u_y, u_z;
	int ncells;
	int nwalls;
	int nnodes;
	int bc[6]; //boundary type
	int flux_bc_flag,cell_bc_flag; //cell_bc_flag is 1 if all the boundaries are updated imposing cell averages
						 //flux_bc_flag is 1 if all the boundaries are updated with numerical fluxes
	int periodicX,periodicY,periodicZ;
	t_cell *cell;
	t_wall *wall;
	t_node *node;

	double lambda_max;
	double tke;
	double mass,energy,mass0,energy0;

	t_sim *sim;

};

struct t_sim_{
	double dt,t,CFL;  //dynamic variables
	double tf, tVolc; //static variables
	int rk_steps; //number of Runge-Kutta steps (1 if 1-st order, 3 otherwise)
	int order; //order of accuracy
	int nvar;  //number of variables

};


struct t_solid_{
	int nsolid;
	char * filename[50];
	t_stl *stl;

};

struct t_stl_{
	int ntri;
	int nver;
	char name[256];
	double Xmin[3],Xmax[3];
	int imin[3],imax[3];
	t_triangle *triangle;

};

struct t_triangle_{
	int outside; // 1 if at least one node is outside of the domain
	double nr[3],absnr;
	double p1[3],p2[3],p3[3];
      int imin[3],imax[3];
};
#endif