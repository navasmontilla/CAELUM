/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2018-2019 The authors.

License type: Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Spain (CC BY-NC-ND 3.0 ES https://creativecommons.org/licenses/by-nc-nd/3.0/es/deed.en) under the following terms:

- Attribution — You must give appropriate credit and provide a link to the license.
- NonCommercial — You may not use the material for commercial purposes.
- NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material unless explicit permission of the authors is provided.

Disclaimer: This software is distributed for research and/or academic purposes, WITHOUT ANY WARRANTY. In no event shall the authors be liable for any claim, damages or other liability, arising from, out of or in connection with the software or the use or other dealings in this software.

File:
  - structures.h

Content:
  -This header file contains the structures designed for the model
  and used in the rest of the code.

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
	double *U; //this is the array of conserved variables. When using Euler: rho, rhou, rhov, E, rhophi
	double *U_aux;
	double *Ue; //equilibrium state
	double *S;
      double *S_corr;
	double pres,prese,u_int; 
	double dx,dy,dz;
	double xc,yc,zc;
	
	int n1,n2,n3,n4,n5,n6,n7,n8; //ID's of the nodes of the cell
	int w1_id,w2_id,w3_id,w4_id,w5_id,w6_id; //ID's of the walls of the cell
	t_wall *w1, *w2, *w3, *w4, *w5, *w6;

      int type; //0=solid, 1=nomal
      int ghost; //0=no, 1=yes
	double xim,yim,zim;
      int ni[8]; //id: for each ghost cell, neighbors of image points.
	double li[8],distabs;
      int distsolx,distsoly,distsolz; //distance to solid cells in cell units
      int solid_id,triangle_id;
      int out; //auxiliary variable for assing_cell_type() . It indicates that a cell is outside of a surface
      t_triangle *tri;
      int st_sizeX, st_sizeY, st_sizeZ;	//stencil size
	int stX[9], stY[9], stZ[9];		//id's of the cells in the X and Y stencil.

};


struct t_wall_{
	int id;
	int stencil; //this could be 1, 3, 5 or 7, depending on the method stencil
	double *UL, *UR; //array of reconstructed values on the left and right hand side of the wall, coming from (WENO/TENO) reconstruction
	double *fR_star,*fL_star;
      double *ULe, *URe; //array of reconstructed values on the left and right hand side of the wall, coming from (WENO/TENO) reconstruction, for the EQUILIBRIUM
	double pRe,pLe; //equilibrium pressures
	int cellR_id, cellL_id; //id of the right and left cell
	t_cell *cellR, *cellL; //pointers to the left and right hand cells of the wall
	double nx, ny, nz;
	double z; //height in z direction
      int wtype, boundId; //wtype: 1 for inner walls, 3 for transmissive boundary walls and 4 for solid walls
      double vel;
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