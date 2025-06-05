/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - preproc.c

Content:
  -This file contains all the pre-processing utilities
  
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>

#include "definitions.h"
#include "structures.h"
#include "preproc.h"
#include "closures.h"
#include "mathutils.h"


int create_mesh(t_mesh *mesh, t_sim *sim){
	int l,m,n,k,aux,k2d;
	int xcells,ycells,zcells;
	t_cell *cell;
	t_wall *wall;
	t_node *node;

	//int semiSt; //auxiliar variable to recalculate boundary cells stencil
	mesh->tke=0.0;


	mesh->sim=sim;

	//Cells
	xcells=mesh->xcells;
	ycells=mesh->ycells;
      zcells=mesh->zcells;
	mesh->ncells=xcells*ycells*zcells;
	mesh->cell=(t_cell*)malloc(mesh->ncells*sizeof(t_cell));
	cell=mesh->cell;

      for(n=0;n<zcells;n++){
            for(m=0;m<ycells;m++){
                  for(l=0;l<xcells;l++){

                        k= xcells*m + l + n*xcells*ycells;
                        cell[k].id=999;
                        cell[k].l=l;
                        cell[k].m=m;
                        cell[k].n=n;


                        cell[k].w1_id=999;
                        cell[k].w2_id=999;
                        cell[k].w3_id=999;
                        cell[k].w4_id=999;
                        cell[k].w5_id=999;
                        cell[k].w6_id=999;

                        cell[k].dx=mesh->dx;
                        cell[k].dy=mesh->dy;
                        cell[k].dz=mesh->dz;

                        cell[k].xc=999;
                        cell[k].yc=999;
                        cell[k].zc=999;

				cell[k].xim = 999;
				cell[k].yim = 999;
				cell[k].zim = 999;

                        cell[k].distabs = 9999999999999.0;

                        cell[k].distsolx = 9999999;
                        cell[k].distsoly = 9999999;
                        cell[k].distsolz = 9999999;

                        cell[k].out = 0;

                        cell[k].n1=999;
                        cell[k].n2=999;
                        cell[k].n3=999;
                        cell[k].n4=999;
                        cell[k].n5=999;
                        cell[k].n6=999;
                        cell[k].n7=999;
                        cell[k].n8=999;
                  }
            }
      }

	//Walls
	mesh->nwalls=3*mesh->ncells+xcells*zcells+ycells*zcells+xcells*ycells;
	mesh->wall=(t_wall*)malloc(mesh->nwalls*sizeof(t_wall));

	wall=mesh->wall;
	for(k=0;k<mesh->nwalls;k++){
		wall[k].id=k;
		wall[k].fR_star=(double*)malloc(sim->nvar*sizeof(double));
		wall[k].fL_star=(double*)malloc(sim->nvar*sizeof(double));
		wall[k].UR=(double*)malloc(sim->nvar*sizeof(double));
		wall[k].UL=(double*)malloc(sim->nvar*sizeof(double));
            wall[k].URe=(double*)malloc(sim->nvar*sizeof(double));
		wall[k].ULe=(double*)malloc(sim->nvar*sizeof(double));
	}

	//Walls amd nodes of the cells
      for(n=0;n<zcells;n++){
            for(m=0;m<ycells;m++){
                  for(l=0;l<xcells;l++){

                        k = l + m*xcells + n*xcells*ycells;
                        k2d=l + m*xcells;
                        cell[k].id=k;
                        cell[k].l=l;
                        cell[k].m=m;


                        cell[k].w1_id=3*(k2d)+m+n*(3*xcells*ycells+xcells+ycells);
                        cell[k].w4_id=cell[k].w1_id+1;
                        cell[k].w5_id=cell[k].w1_id+2;
                        if (l==xcells-1){
                              cell[k].w2_id=cell[k].w1_id+4 - 1;
                        }else{
                              cell[k].w2_id=cell[k].w1_id+4;
                        }
                        if (m==ycells-1){
                              aux=(3*xcells*ycells+xcells+ycells)*(n+1);
                              cell[k].w3_id=aux-xcells+l;
                        }else{
                              cell[k].w3_id=cell[k].w1_id+xcells*3+1;
                        }
                        if (n==zcells-1){
                              cell[k].w6_id=mesh->nwalls-xcells*ycells+l+m*xcells;
                        }else{
                              cell[k].w6_id=cell[k].w5_id+(3*xcells*ycells+xcells+ycells);
                        }

                        cell[k].w1=&(mesh->wall[cell[k].w1_id]);
                        cell[k].w2=&(mesh->wall[cell[k].w2_id]);
                        cell[k].w3=&(mesh->wall[cell[k].w3_id]);
                        cell[k].w4=&(mesh->wall[cell[k].w4_id]);
                        cell[k].w5=&(mesh->wall[cell[k].w5_id]);
                        cell[k].w6=&(mesh->wall[cell[k].w6_id]);

                        cell[k].dx=mesh->dx;
                        cell[k].dy=mesh->dy;
                        cell[k].dz=mesh->dz;

                        cell[k].xc=(l+0.5)*cell[k].dx;
                        cell[k].yc=(m+0.5)*cell[k].dy;
                        cell[k].zc=(n+0.5)*cell[k].dz;

                        cell[k].n1=k2d +m +n*(xcells+1)*(ycells+1);
                        cell[k].n2=cell[k].n1+1;
                        cell[k].n3=cell[k].n2+xcells+1;
                        cell[k].n4=cell[k].n2+xcells;
                        cell[k].n5=cell[k].n1 + (xcells+1)*(ycells+1);
                        cell[k].n6=cell[k].n5+1;
                        cell[k].n7=cell[k].n6+xcells+1;
                        cell[k].n8=cell[k].n6+xcells;

                  }
            }
      }


	//normal vectors of the walls
  	for(k=0;k<mesh->ncells;k++){
		mesh->wall[cell[k].w1_id].nx=0.0;
		mesh->wall[cell[k].w4_id].nx=1.0;
		mesh->wall[cell[k].w1_id].ny=1.0;
		mesh->wall[cell[k].w4_id].ny=0.0;
            mesh->wall[cell[k].w1_id].nz=0.0;
		mesh->wall[cell[k].w4_id].nz=0.0;

		//we know this is redundant,
		//but it's easier and is pre-process
		mesh->wall[cell[k].w3_id].nx=0.0;
		mesh->wall[cell[k].w2_id].nx=1.0;
		mesh->wall[cell[k].w3_id].ny=1.0;
		mesh->wall[cell[k].w2_id].ny=0.0;
            mesh->wall[cell[k].w3_id].nz=0.0;
		mesh->wall[cell[k].w2_id].nz=0.0;

            mesh->wall[cell[k].w5_id].nx=0.0;
		mesh->wall[cell[k].w6_id].nx=0.0;
		mesh->wall[cell[k].w5_id].ny=0.0;
		mesh->wall[cell[k].w6_id].ny=0.0;
            mesh->wall[cell[k].w5_id].nz=1.0;
		mesh->wall[cell[k].w6_id].nz=1.0;

	}



	//Nodes
	mesh->nnodes=(mesh->xcells+1)*(mesh->ycells+1)*(mesh->zcells+1);
	mesh->node=(t_node*)malloc(mesh->nnodes*sizeof(t_node));

	node=mesh->node;
	//Loop for nodes
      for(n=0;n<zcells+1;n++){
            for(m=0;m<ycells+1;m++){
                  for(l=0;l<xcells+1;l++){
                        k=(xcells+1)*(ycells+1)*n+(xcells+1)*m+l;
                        node[k].id=k;
                        node[k].x=l*mesh->dx;
                        node[k].y=m*mesh->dy;
                        node[k].z=n*mesh->dz;
                  }
            }
      }

	//Boundary condition flags
	//This is for BC imposition
	if(mesh->bc[0]==99 && mesh->bc[1]==99 && mesh->bc[2]==99 && mesh->bc[3]==99){ //this has been temporary deactivated
		//all the boundaries are updated imposing cell averages
		mesh->cell_bc_flag=1;
 	}else{
		mesh->cell_bc_flag=0;
	}
	if(mesh->bc[0]!=1 || mesh->bc[1]!=1 || mesh->bc[2]!=1 || mesh->bc[3]!=1 || mesh->bc[4]!=1 || mesh->bc[5]!=1){
		//all the boundaries are updated with fluxes
		mesh->flux_bc_flag=1;
 	}else{
		mesh->flux_bc_flag=0;
	}

	//This defines boundaries for interpolation
	if(mesh->bc[1]==1 && mesh->bc[3]==1){
		mesh->periodicX=1;
	}else{
		mesh->periodicX=0;
	}
	if(mesh->bc[0]==1 && mesh->bc[2]==1){
		mesh->periodicY=1;
	}else{
		mesh->periodicY=0;
	}
      if(mesh->bc[4]==1 && mesh->bc[5]==1){
		mesh->periodicZ=1;
	}else{
		mesh->periodicZ=0;
	}


	//Assigment of wall's neighbour cells
      for(n=0;n<zcells;n++){
            for(m=0;m<ycells;m++){
                  for(l=0;l<xcells;l++){

                        //Interior walls
                        k = l + m*xcells + n*xcells*ycells;

                        cell[k].w1->cellR_id=cell[k].id;
                        cell[k].w4->cellR_id=cell[k].id;
                        cell[k].w5->cellR_id=cell[k].id;

                        cell[k].w2->cellL_id=cell[k].id;
                        cell[k].w3->cellL_id=cell[k].id;
                        cell[k].w6->cellL_id=cell[k].id;

                        cell[k].w1->cellR=&(cell[k]);
                        cell[k].w4->cellR=&(cell[k]);
                        cell[k].w5->cellR=&(cell[k]);
                        cell[k].w2->cellL=&(cell[k]);
                        cell[k].w3->cellL=&(cell[k]);
                        cell[k].w6->cellL=&(cell[k]);

                        //Boundary walls need to stablish their external cells
                        if(m==0){
                              //If we are at the first row position (m=0)
                              //the left cell of the boundary walls are the last row cells
                              //assuming periodic boundaries (if not, the fluxes will be rewritten later on
                              cell[k].w1->cellL_id=cell[k+(ycells-1)*xcells].id;  //xcells*(ycells-1)+l+n*zcell
                              cell[k].w1->cellL=&(cell[k+(ycells-1)*xcells]);
                        }

                        if(m==ycells-1){
                              //If we are at the last row position (m=ycells-1)
                              //the right cell of the boundary walls are the first row cells
                              //assuming periodic boundaries (if not, the fluxes will be rewritten later on
                              cell[k].w3->cellR_id=cell[k-m*xcells].id;
                              cell[k].w3->cellR=&(cell[k-m*xcells]);
                        }

                        if(l==0){
                              cell[k].w4->cellL_id=cell[k+(xcells-1)].id;
                              cell[k].w4->cellL=&(cell[k+(xcells-1)]);
                        }
                        if(l==xcells-1){
                              cell[k].w2->cellR_id=cell[k-(xcells-1)].id;
                              cell[k].w2->cellR=&(cell[k-(xcells-1)]);
                        }

                        if(n==0){
                              cell[k].w5->cellL_id=cell[k+(zcells-1)*xcells*ycells].id;
                              cell[k].w5->cellL=&(cell[k+(zcells-1)*xcells*ycells]);
                        }
                        if(n==zcells-1){
                              cell[k].w6->cellR_id=cell[k-(zcells-1)*xcells*ycells].id;
                              cell[k].w6->cellR=&(cell[k-(zcells-1)*xcells*ycells]);
                        }

                  }
            }
      }
	//Allocation of arrays of variables in cells and walls


	for(k=0;k<mesh->ncells;k++){
		mesh->cell[k].U=    (double*)malloc(sim->nvar*sizeof(double));
		mesh->cell[k].U_aux=(double*)malloc(sim->nvar*sizeof(double));
		mesh->cell[k].Ue   =(double*)malloc(sim->nvar*sizeof(double));
		mesh->cell[k].S=    (double*)malloc(sim->nvar*sizeof(double));
            mesh->cell[k].S_corr=    (double*)malloc(sim->nvar*sizeof(double));
	}

	for(k=0;k<mesh->nwalls;k++){
		mesh->wall[k].UR=(double*)malloc(sim->nvar*sizeof(double));
		mesh->wall[k].UL=(double*)malloc(sim->nvar*sizeof(double));
	}

	cell=mesh->cell;
	for(n=0;n<mesh->ncells;n++){
		for(k=0;k<sim->nvar;k++){
			cell[n].S[k]=0.0 ;
                  cell[n].S_corr[k]=0.0 ;
		}
	}

	for(n=0;n<mesh->ncells;n++){
		cell[n].w1->z=cell[n].zc;
		cell[n].w4->z=cell[n].zc;
		cell[n].w2->z=cell[n].zc;
		cell[n].w3->z=cell[n].zc;
		cell[n].w5->z=cell[n].zc-0.5*cell[n].dz;
		cell[n].w6->z=cell[n].zc+0.5*cell[n].dz;
	}


	printf("%s Memory has been allocated and mesh connectivity has been defined \n",OK);

	return 1;
}


int read_initial(t_mesh *mesh, t_sim *sim, const char *folder_path){
	int m,k,l,n,ct;
	FILE *fp;
	char fname[1024];
      double u;
#if EQUATION_SYSTEM == 2
      double p,v,w,rho,phi,gamma;
      #if ST!=0
      FILE *fpe;
      #endif
#endif
	
	t_cell *cell;
	cell=mesh->cell;
	ct=0;
	
	snprintf(fname, sizeof(fname), "%s/initial.out", folder_path);
	fp = fopen(fname,"r");

#if EQUATION_SYSTEM == 2
	#if ST!=0
      
      snprintf(fname, sizeof(fname), "%s/equilibrium.out", folder_path);
	fpe = fopen(fname,"r");
      
	if (fpe != NULL){
		// Skip the first two lines
		if (fscanf(fpe, "%*[^\n]\n") != 0) {
		  printf("Warning: Failed to skip the first line.\n");
		}
		if (fscanf(fpe, "%*[^\n]\n") != 0) {
		  printf("Warning: Failed to skip the second line.\n");
		}

		for(l=0;l<mesh->xcells;l++){
			for(m=0;m<mesh->ycells;m++){
				for(n=0;n<mesh->zcells;n++){
				k = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
				if (fscanf(fpe, "%*f %*f %*f %le %le %le %le %le %le", &u, &v, &w, &rho, &p, &phi) != 6) {
				printf("%s Error: Failed to read data equilibrium data \n",WAR);
				getchar();
				}
				#if MULTICOMPONENT
					#if MULTI_TYPE==1
					gamma=phi;
					#else
					gamma=1.0+1.0/phi;
					#endif
				#else
					gamma=_gamma_;
				#endif
				cell[k].Ue[0]=rho;
				cell[k].Ue[1]=u*cell[k].U[0];
				cell[k].Ue[2]=v*cell[k].U[0];
				cell[k].Ue[3]=w*cell[k].U[0];
				cell[k].Ue[4]=energy_from_pressure(gamma,p,u,v,w,rho,cell[k].zc);
				cell[k].Ue[5]=phi*rho;
				cell[k].prese=p;
				}
			}
		}

		fclose(fpe);

		printf("%s equilibrium.out file has been read \n",OK);
	}else{
		printf("%s File equilibrium.out not found. Initial and equilibrium data is set in update_initial() \n",WAR);
		ct=1;
	}
	#endif
	
	if (fp != NULL){
		// Skip the first two lines
		if (fscanf(fp, "%*[^\n]\n") != 0) {
		  printf("Warning: Failed to skip the first line.\n");
		}
		if (fscanf(fp, "%*[^\n]\n") != 0) {
		  printf("Warning: Failed to skip the second line.\n");
		}

		for(l=0;l<mesh->xcells;l++){
			for(m=0;m<mesh->ycells;m++){
				for(n=0;n<mesh->zcells;n++){
				k = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
				if (fscanf(fp, "%*f %*f %*f %le %le %le %le %le %le", &u, &v, &w, &rho, &p, &phi) != 6) {
				printf("%s Error: Failed to read data initial data \n",WAR);
				getchar();
				}
				#if MULTICOMPONENT
					#if MULTI_TYPE==1
					gamma=phi;
					#else
					gamma=1.0+1.0/phi;
					#endif
				#else
					gamma=_gamma_;
				#endif
				cell[k].U[0]=rho;
				cell[k].U[1]=u*cell[k].U[0];
				cell[k].U[2]=v*cell[k].U[0];
				cell[k].U[3]=w*cell[k].U[0];
				cell[k].U[4]=energy_from_pressure(gamma,p,u,v,w,rho,cell[k].zc);
				cell[k].U[5]=phi*rho;
				}
			}
		}

		fclose(fp);

		printf("%s initial.out file has been read \n",OK);

	}else{

	printf("%s File initial.out not found. Initial and equilibrium data is set in update_initial() \n",WAR);
	ct=1;
	
	}
	
#else 
	if (fp != NULL){
		// Skip the first two lines
		if (fscanf(fp, "%*[^\n]\n") != 0) {
		  printf("Warning: Failed to skip the first line.\n");
		}
		if (fscanf(fp, "%*[^\n]\n") != 0) {
		  printf("Warning: Failed to skip the second line.\n");
		}

		for(l=0;l<mesh->xcells;l++){
			for(m=0;m<mesh->ycells;m++){
				for(n=0;n<mesh->zcells;n++){
				k = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
				if (fscanf(fp, "%*f %*f %*f %le ", &u) != 1) {
				printf("%s Error: Failed to read data initial data \n",WAR);
				getchar();
				}
				cell[k].U[0]=u;
				}
			}
		}

		fclose(fp);

		printf("%s initial.out file has been read \n",OK);

	}else{

	printf("%s File initial.out not found. Initial and equilibrium data is set in update_initial() \n",WAR);
	ct=1;
	
	}

#endif
	
	return ct;
	
}



int update_initial(t_mesh *mesh, t_sim *sim, const char *folder_path){
	int k,ct;
      double xc;
#if EQUATION_SYSTEM == 2 
      int m;
      double zc,d1,d2,rc,aux1,aux2;
      double p,u,v,w,rho,phi,gamma,tt,p0,tt0,rho0;
#else
      double yc,r;
#endif

	t_cell *cell;
	cell=mesh->cell;

#if READ_INITIAL
	ct=read_initial(mesh,sim,folder_path);
	if(ct==1){
#else
	printf("%s Read initial data from file case/initial.out is disabled \n",WAR);
#endif

#if EQUATION_SYSTEM == 2 

      for(k=0;k<mesh->ncells;k++){


            if(cell[k].type==1){

                  u=0.0;
                  v=0.0;
                  w=0.0;
                  phi=0.0;


                  tt0=300;
                  p0=_p0_;
                  rho0=p0/(_R_*tt0);

                  tt=tt0;
                  
                  aux2=(_gamma_-1.0)/_gamma_*_g_/(_R_*tt);


                  p=   p0*  pow((1.0-aux2*cell[k].zc),_gamma_/(_gamma_-1.0));  //pointwise
                  rho= rho0*pow((1.0-aux2*cell[k].zc),1.0/(_gamma_-1.0));


			gamma=_gamma_;
			cell[k].Ue[0]=rho;
			cell[k].Ue[1]=u*cell[k].Ue[0];
			cell[k].Ue[2]=v*cell[k].Ue[0];
			cell[k].Ue[3]=w*cell[k].Ue[0];
			cell[k].Ue[4]=energy_from_pressure(gamma,p,u,v,w,rho,cell[k].zc);
			cell[k].Ue[5]=phi;

                  cell[k].prese=p;

			for(m=0;m<sim->nvar;m++){
				cell[k].U[m]=cell[k].Ue[m];
			}


                  xc=10000;
                  zc=2000;
                  d1=sqrt((cell[k].xc-xc)*(cell[k].xc-xc)+(cell[k].zc-zc)*(cell[k].zc-zc));

                  xc=10000;
                  zc=8000;
                  d2=sqrt((cell[k].xc-xc)*(cell[k].xc-xc)+(cell[k].zc-zc)*(cell[k].zc-zc));

                  rc=1000;

                  aux1=20.0*(MAX(rc-d1/2.0,0.0)+MIN(d2/2.0-rc,0.0))/1000;


                  tt=tt0+aux1;
                  aux2=(_gamma_-1.0)/_gamma_*_g_/(_R_*tt0);
                  p=   p0*  pow((1.0-aux2*cell[k].zc),_gamma_/(_gamma_-1.0));  //pointwise
                  rho= p0/(_R_*tt)*pow((1.0-aux2*cell[k].zc),1.0/(_gamma_-1.0));
			u=0.0;
                  v=0.0;
                  w=0.0;
                  phi=0.0;

			gamma=_gamma_;
			cell[k].U[0]=rho;
			cell[k].U[1]=u*cell[k].U[0];
			cell[k].U[2]=v*cell[k].U[0];
			cell[k].U[3]=w*cell[k].U[0];
			cell[k].U[4]=energy_from_pressure(gamma,p,u,v,w,rho,cell[k].zc);
			cell[k].U[5]=phi;


            }else{

                  cell[k].Ue[0]=-1.0;
                  cell[k].Ue[1]=0.0;
                  cell[k].Ue[2]=0.0;
                  cell[k].Ue[3]=0.0;
                  cell[k].Ue[4]=0.0;
                  cell[k].Ue[5]=0.0;

			for(m=0;m<sim->nvar;m++){
				cell[k].U[m]=cell[k].Ue[m];
			}


            }

      }

#if READ_INITIAL
	}
#endif

#endif



#if EQUATION_SYSTEM == 1
for(k=0;k<mesh->ncells;k++){

      xc=5.0;
      yc=5.0;
      //zc=5.0;
      //r=sqrt((cell[k].xc-xc)*(cell[k].xc-xc)+(cell[k].yc-yc)*(cell[k].yc-yc)+(cell[k].zc-zc)*(cell[k].zc-zc));
	r=sqrt((cell[k].xc-xc)*(cell[k].xc-xc)+(cell[k].yc-yc)*(cell[k].yc-yc));

      if (r<0.5) {
            cell[k].U[0]=2.0;
      }else{
            cell[k].U[0]=1.0;
      }

}

#if READ_INITIAL
	}
#endif

#endif

#if EQUATION_SYSTEM == 0 
for(k=0;k<mesh->ncells;k++){

      xc=5.0;
      yc=5.0;
      //zc=5.0;
      //r=sqrt((cell[k].xc-xc)*(cell[k].xc-xc)+(cell[k].yc-yc)*(cell[k].yc-yc)+(cell[k].zc-zc)*(cell[k].zc-zc));
	r=sqrt((cell[k].xc-xc)*(cell[k].xc-xc)+(cell[k].yc-yc)*(cell[k].yc-yc));

      if (r<0.5) {
            cell[k].U[0]=2.0;
      }else{
            cell[k].U[0]=1.0;
      }

}

#if READ_INITIAL
	}
#endif

#endif


	return 1;
}



int assign_wall_type(t_mesh *mesh){
	t_wall *wall;
	int m,l,n,k;

    for(n=0;n<mesh->nwalls;n++){
            wall=&(mesh->wall[n]);
            wall->wtype=1;          //by default: 1= normal RP wall
            wall->boundId=999;    //999 when the wall is not at any boundary. Otherwise: 1, 2, 3, 4, 5, 6.

            #if ALLOW_SOLIDS==1


            if(wall->nx>TOL4){
                  if(wall->cellL->type==0 && wall->cellR->type!=0){
                        wall->wtype=4;
                        wall->boundId=4;
                  }else if(wall->cellR->type==0 && wall->cellL->type!=0){
                        wall->wtype=4;
                        wall->boundId=2;
                  }else if(wall->cellL->type==0 && wall->cellR->type==0){
                        wall->wtype=0;
                  }
            }else if(wall->ny>TOL4){
                  if(wall->cellL->type==0 && wall->cellR->type!=0){
                        wall->wtype=4;
                        wall->boundId=1;
                  }else if(wall->cellR->type==0 && wall->cellL->type!=0){
                        wall->wtype=4;
                        wall->boundId=3;
                  }else if(wall->cellL->type==0 && wall->cellR->type==0){
                        wall->wtype=0;
                  }
            }else{
                  if(wall->cellL->type==0 && wall->cellR->type!=0){
                        wall->wtype=4;
                        wall->boundId=5;
                  }else if(wall->cellR->type==0 && wall->cellL->type!=0){
                        wall->wtype=4;
                        wall->boundId=6;
                  }else if(wall->cellL->type==0 && wall->cellR->type==0){
                        wall->wtype=0;
                  }
            }



            #endif

    }



      m=0;
      for(l=0;l<mesh->xcells;l++){
            for(n=0;n<mesh->zcells;n++){
                  //k=mesh->xcells*m+l;
                  k= mesh->xcells*m + l + n*mesh->xcells*mesh->ycells;
                  wall=&(mesh->wall[mesh->cell[k].w1_id]);
                  if(mesh->cell[k].type!=0){
                        wall->wtype=mesh->bc[0];
                        wall->boundId=1;
                  }else{
                        wall->wtype=0;
                  }

            }
      }

      m=mesh->ycells-1;
      for(l=0;l<mesh->xcells;l++){
            for(n=0;n<mesh->zcells;n++){
                  k= mesh->xcells*m + l + n*mesh->xcells*mesh->ycells;
                  wall=&(mesh->wall[mesh->cell[k].w3_id]);
                  if(mesh->cell[k].type!=0){
                        wall->wtype=mesh->bc[2];
                        wall->boundId=3;
                  }else{
                        wall->wtype=0;
                  }
            }
      }


      l=mesh->xcells-1;
      for(m=0;m<mesh->ycells;m++){
            for(n=0;n<mesh->zcells;n++){
                  k= mesh->xcells*m + l + n*mesh->xcells*mesh->ycells;
                  wall=&(mesh->wall[mesh->cell[k].w2_id]);
                  if(mesh->cell[k].type!=0){
                        wall->wtype=mesh->bc[1];
                        wall->boundId=2;
                  }else{
                        wall->wtype=0;
                  }
            }
      }


      l=0;
      for(m=0;m<mesh->ycells;m++){
            for(n=0;n<mesh->zcells;n++){
                  k= mesh->xcells*m + l + n*mesh->xcells*mesh->ycells;
                  wall=&(mesh->wall[mesh->cell[k].w4_id]);
                  if(mesh->cell[k].type!=0){
                        wall->wtype=mesh->bc[3];
                        wall->boundId=4;
                  }else{
                        wall->wtype=0;
                  }
            }
      }

      n=0;
      for(m=0;m<mesh->ycells;m++){
            for(l=0;l<mesh->xcells;l++){
                  k= mesh->xcells*m + l + n*mesh->xcells*mesh->ycells;
                  wall=&(mesh->wall[mesh->cell[k].w5_id]);
                  if(mesh->cell[k].type!=0){
                        wall->wtype=mesh->bc[4];
                        wall->boundId=5;
                  }else{
                        wall->wtype=0;
                  }
            }
      }

      n=mesh->zcells-1;;
      for(m=0;m<mesh->ycells;m++){
            for(l=0;l<mesh->xcells;l++){
                  k= mesh->xcells*m + l + n*mesh->xcells*mesh->ycells;
                  wall=&(mesh->wall[mesh->cell[k].w6_id]);
                  if(mesh->cell[k].type!=0){
                        wall->wtype=mesh->bc[5];
                        wall->boundId=6;
                  }else{
                        wall->wtype=0;
                  }
            }
      }


	return 1;
}

int assign_cell_type(t_mesh *mesh,t_solid *solids){ // Define ghost and solid cells
	t_cell *cell;	
	int k;
#if ALLOW_SOLIDS==1
	int m,l,n,i,j,q,i1,i2,i3,na,df,df0,ct;
      double aux1,aux2;
      double proj,dist,s1,s2,s3;
      double dif[3],xc[3],v1[3],v2[3],v3[3],vp1[3],vp2[3],vp3[3],dc1[3],dc2[3],dc3[3];
      double dp;
      int solx,soly,solz;
      t_triangle *triangle;
#endif

      cell=mesh->cell;
      for(k=0;k<mesh->ncells;k++){
            cell[k].type=1;          //by default 1.    1= computed cell, 0= solid cell (not computed cell)
            cell[k].ghost=0;
      }

#if ALLOW_SOLIDS==1
      /*
      //solid cells are assigned using simple formulas. In the future, "find-point-inside" algorithms will be used.
      */

      if(solids->nsolid<1){
            printf("%s In function assign_cell_type() no solids are considered\n",WAR);
      }else{

            for(l=0;l<solids->nsolid;l++){ // loop for each solid
                  triangle=solids->stl[l].triangle;
                  //printf("Solid %d: boundary cells are being computed...\n\n",l);
                  for(m=0;m<solids->stl[l].ntri;m++){ // loop for each triangle
                        //if(triangle[m].outside==0){
                              for(q=0;q<3;q++){
                                    v1[q]=triangle[m].p2[q]-triangle[m].p1[q]; // the vectors between nodes are defined
                                    v2[q]=triangle[m].p3[q]-triangle[m].p2[q];
                                    v3[q]=triangle[m].p1[q]-triangle[m].p3[q];
                              }
                              dp=MAX(mesh->dx,mesh->dy);
                              dp=MAX(dp,mesh->dz);
                              dp=_stol_*dp ; //2.0*fabs(cell[n].dx*triangle[m].nr[0]+cell[n].dy*triangle[m].nr[1]+cell[n].dz*triangle[m].nr[2])/triangle[m].absnr;
                              //dp=_stol_*mesh->dx;

                              //loop over cells inside the bounding box of the triangle
                              for(i=triangle[m].imin[0];i<=triangle[m].imax[0];i++){
                                    for(j=triangle[m].imin[1];j<=triangle[m].imax[1];j++){
                                          for(k=triangle[m].imin[2];k<=triangle[m].imax[2];k++){
                                                //printf("triangle %d, cell (%d %d %d) \n",m,i,j,k);
                                                n= mesh->xcells*j + i + k*mesh->xcells*mesh->ycells;
                                                dif[0]=cell[n].xc-triangle[m].p1[0]; //vector from cell center to node P1 (x-component)
                                                dif[1]=cell[n].yc-triangle[m].p1[1]; //vector from cell center to node P1 (y-component)
                                                dif[2]=cell[n].zc-triangle[m].p1[2]; //vector from cell center to node P1 (z-component)
                                                proj=dif[0]*triangle[m].nr[0]+dif[1]*triangle[m].nr[1]+dif[2]*triangle[m].nr[2]; //vector from cell center to triangle P1, projected onto normal direction
                                                dist=proj/triangle[m].absnr;
                                                if (proj>0 && cell[n].ghost>0 && fabs(dist)<cell[n].distabs){ // if a cell has been tagged as "ghost" and is now found outside of other triangle, it is reverted as "fluid cell"
                                                      cell[n].ghost=0;
                                                }
                                                if (proj<=0&&cell[n].out<1) { //when proj<0 and the cell is not already outside of other triangle, the cell center is inside the solid (below the surface)
                                                      if (fabs(dist)<dp) {
                                                            xc[0]=cell[n].xc - dist*triangle[m].nr[0]; //Formula: Xc = Xcell - proj * n / abs(n) allows to compute the intersection point on the triangle, xc
                                                            xc[1]=cell[n].yc - dist*triangle[m].nr[1];
                                                            xc[2]=cell[n].zc - dist*triangle[m].nr[2];

                                                            // Below, the point-in-triangle algorithm (3D version) is used to detect if the intersection point, xc, is inside the triangle
                                                            for(q=0;q<3;q++){
                                                                  dc1[q]=xc[q]-triangle[m].p1[q];
                                                                  dc2[q]=xc[q]-triangle[m].p2[q];
                                                                  dc3[q]=xc[q]-triangle[m].p3[q];
                                                            }

                                                            vector_product(v1,dc1,vp1); // the point-in-triangle algorithm  is based on cross products to determine whether point is inside or not
                                                            vector_product(v2,dc2,vp2);
                                                            vector_product(v3,dc3,vp3);

                                                            s1=dot_product(vp1,vp2);
                                                            s2=dot_product(vp2,vp3);
                                                            s3=dot_product(vp3,vp1);

                                                            if (s1>0 && s2>0 &&  s3>0) { // when all are positive (same vector product direction), the point xc is inside the triangle
                                                                  //cell[n].type=0;
                                                                  if(fabs(dist)<cell[n].distabs){  // the condition before may happen for many triangles, so we keep the closest to the interface
                                                                        cell[n].ghost=1;
                                                                        cell[n].solid_id=l;
                                                                        cell[n].triangle_id=m;
                                                                        cell[n].tri=&(triangle[m]);
                                                                        if(fabs(dist)<TOL14){
                                                                              dist=-TOL14; //to ensure negativity of dist, that might be zero and therefore +-0
                                                                        }
                                                                        cell[n].distabs=fabs(dist);
                                                                        //printf("triangle %14.14e \n",dist);
                                                                        cell[n].xim = xc[0] - dist*triangle[m].nr[0];
                                                                        cell[n].yim = xc[1] - dist*triangle[m].nr[1];
                                                                        cell[n].zim = xc[2] - dist*triangle[m].nr[2];
                                                                  }
                                                            }
                                                      }

                                                }else{
                                                            cell[n].out=1;

                                                }

                                          }
                                    }
                              }
                              //getchar();

                        //}

                  }
            }

            //printf("%s Ghost cells have been identified \n\n",OK);


            for(l=0;l<solids->nsolid;l++){ //loop over cells inside the solid bounding box
                  for(i=solids->stl[l].imin[0];i<=solids->stl[l].imax[0];i++){
                              for(j=solids->stl[l].imin[1];j<=solids->stl[l].imax[1];j++){
                                    for(k=solids->stl[l].imin[2];k<=solids->stl[l].imax[2];k++){

                                          n= mesh->xcells*j + i + k*mesh->xcells*mesh->ycells;
                                          if(cell[n].ghost!=1){

                                                solx=0;
                                                soly=0;
                                                solz=0;

                                                df0=mesh->xcells;
                                                q=-1;
                                                ct=0; // we define a counter and only if ct is greater or equal to 2, solids are defined, to prevent from individual ghost cells generating solid lines
                                                for (i1=0;i1<mesh->xcells;i1++){ //loop over lines to determine the closest interface to a point, in x-direction
                                                      na=mesh->xcells*j + i1 + k*mesh->xcells*mesh->ycells;
                                                      if(cell[na].ghost==1){
                                                            ct+=1; //counter to check the number of ghost cells in a line
                                                            df=abs(i1-i); //distance from i1 to i
                                                            if(df<=df0){
                                                                  df0=df;
                                                                  q=na; //this is the index of the nearest ghost cell to cell n
                                                            }
                                                      }
                                                }
                                                if(q>-1&&ct>1){
                                                      aux1=cell[q].xc-cell[n].xc;
                                                      aux2=aux1*cell[q].tri->nr[0];
                                                      if(aux2>0.0){ //when the vector points out the surface, the cell is inside
                                                            solx=1;//cell[n].type=0;
                                                      }
                                                }

                                                //if(cell[n].xc>0.75&&cell[n].type==0){
                                                //      printf("cell %d: q is %d, located at x=%lf \n",n,q,cell[q].xc);
                                                //      printf("ny=%lf \n",cell[q].tri->nr[0]);
                                                //      getchar();
                                                //}

                                                //same than above, for "y"
                                                if(cell[n].type!=0){
                                                df0=mesh->ycells;
                                                q=-1;
                                                ct=0;
                                                for (i2=0;i2<mesh->ycells;i2++){
                                                      na=mesh->xcells*i2 + i + k*mesh->xcells*mesh->ycells;
                                                      if(cell[na].ghost==1){
                                                            ct+=1;
                                                            df=abs(i2-j);
                                                            if(df<=df0){
                                                                  df0=df;
                                                                  q=na;
                                                            }
                                                      }
                                                }
                                                if(q>-1&&ct>1){
                                                      aux1=cell[q].yc-cell[n].yc;
                                                      aux2=aux1*cell[q].tri->nr[1];
                                                      if(aux2>0.0){
                                                            soly=1;//cell[n].type=0;
                                                      }
                                                }
                                                }

                                                //same than above, for "z"
                                                if(cell[n].type!=0){
                                                df0=mesh->zcells;
                                                q=-1;
                                                ct=0;
                                                for (i3=0;i3<mesh->zcells;i3++){
                                                      na=mesh->xcells*j + i + i3*mesh->xcells*mesh->ycells;
                                                      if(cell[na].ghost==1){
                                                            ct+=1;
                                                            df=abs(i3-k);
                                                            if(df<=df0){
                                                                  df0=df;
                                                                  q=na;
                                                            }
                                                      }
                                                }
                                                if(q>-1&&ct>1){
                                                      aux1=cell[q].zc-cell[n].zc;
                                                      aux2=aux1*cell[q].tri->nr[2];
                                                      if(aux2>0.0){
                                                            solz=1;//cell[n].type=0;
                                                      }
                                                }
                                                }

                                                if(solx==1&&soly==1){ //when intersections are detected in X and Y ray tracing, then the cell is set as solid.
                                                      cell[n].type=0;
                                                }


                                          }

                                    }

                              }

                  }


            }

            //looking for orphan cells
            for(l=0;l<solids->nsolid;l++){ //loop over cells inside the solid bounding box
                  for(i=solids->stl[l].imin[0];i<=solids->stl[l].imax[0];i++){
                              for(j=solids->stl[l].imin[1];j<=solids->stl[l].imax[1];j++){
                                    for(k=solids->stl[l].imin[2];k<=solids->stl[l].imax[2];k++){

                                          n= mesh->xcells*j + i + k*mesh->xcells*mesh->ycells;
                                          if(cell[n].ghost!=1){

                                          ct=0;
                                          if(cell[n].type==0){
                                                if(cell[n].w1->cellL->type==0){
                                                      ct+=1;
                                                }
                                                if(cell[n].w2->cellR->type==0){
                                                      ct+=1;
                                                }
                                                if(cell[n].w3->cellR->type==0){
                                                      ct+=1;
                                                }
                                                if(cell[n].w4->cellL->type==0){
                                                      ct+=1;
                                                }
                                                if(cell[n].w5->cellL->type==0){
                                                      ct+=1;
                                                }
                                                if(cell[n].w6->cellR->type==0){
                                                      ct+=1;
                                                }
                                                //if(ct<2||cell[n].out==1){ //when a solid cell has only ONE solid neighbor OR it is ouside the surface, it is re-converted to fluid cell.
                                                if(ct<2){ //when a solid cell has only ONE solid neighbor, it is re-converted to fluid cell.
                                                      cell[n].type=1;
                                                }

                                          }

                                          }


                                    }
                        }
                  }
            }

            //printf("%s Solid cells have been identified \n\n",OK);


      // distance from all cells to the closest solid cell, in the cartesian directions, is computed
      for(k=0;k<mesh->ncells;k++){
            if(cell[k].type!=0){
            m=cell[k].m;
            n=cell[k].n;
            for (l=0;l<mesh->xcells;l++){
                  na = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
                  if(cell[na].type==0){
                        cell[k].distsolx=MIN(cell[k].distsolx,abs(cell[k].l-l))  ;
                  }
            }
            l=cell[k].l;
            n=cell[k].n;
            for (m=0;m<mesh->ycells;m++){
                  na = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
                  if(cell[na].type==0){
                        cell[k].distsoly=MIN(cell[k].distsoly,abs(cell[k].m-m)) ;
                  }
            }
            l=cell[k].l;
            m=cell[k].m;
            for (n=0;n<mesh->zcells;n++){
                  na = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
                  if(cell[na].type==0){
                        cell[k].distsolz=MIN(cell[k].distsolz,abs(cell[k].n-n)) ;
                  }
            }
            //if(cell[k].distsoly<100){
            //      printf("celda %d, dist: %d \n",k,cell[k].distsoly);
            //}

            }
      }


      }

	if(solids->nsolid>0){
      printf("%s Ghost and solid cells have been identified \n",OK);}
	
#endif
	return 1;
}

int update_stencils(t_mesh *mesh,t_sim *sim){
	int l,m,n,k,p,k2d;
	int xcells,ycells,zcells;
	t_cell *cell;

	int semiSt; //auxiliar variable to recalculate boundary cells stencil

	//Cells
	xcells=mesh->xcells;
	ycells=mesh->ycells;
      zcells=mesh->zcells;
	cell=mesh->cell;


      //Set cell stencils
	//Initially all the cells have a stencil of size order
	for(k=0;k<mesh->ncells;k++){
		cell[k].st_sizeX=sim->order;
		cell[k].st_sizeY=sim->order;
            cell[k].st_sizeZ=sim->order;
	}

	//But there are special cases at boundary cells
	semiSt=(sim->order-1)*0.5;
      for(n=0;n<zcells;n++){
            for(m=0;m<ycells;m++){
                  for(l=0;l<xcells;l++){
                        k = l + m*xcells + n*xcells*ycells;
                        k2d=l + m*xcells;
                        if(mesh->periodicX==0){
                              //x setencils
                              if(l<semiSt){
                                    cell[k].st_sizeX=MIN(cell[k].st_sizeX,2*l+1); //poner if elseif else y añadir en todas min(st_size,dist*2+1)
                              }else if(xcells-(l+1)<semiSt){
                                    cell[k].st_sizeX=MIN(cell[k].st_sizeX,2*(xcells-(l+1))+1);
                              }
                              cell[k].st_sizeX=MIN(cell[k].st_sizeX,2*cell[k].distsolx-1);

                        }
                        if(mesh->periodicY==0){
                              //y stencils
                              if(m<semiSt){
                                    cell[k].st_sizeY=MIN(cell[k].st_sizeY,2*m+1);
                              }else if(ycells-(m+1)<semiSt){
                                    cell[k].st_sizeY=MIN(cell[k].st_sizeY,2*(ycells-(m+1))+1);
                              }
                              cell[k].st_sizeY=MIN(cell[k].st_sizeY,2*cell[k].distsoly-1);
                        }
                        if(mesh->periodicZ==0){
                              //y stencils
                              if(n<semiSt){
                                    cell[k].st_sizeZ=MIN(cell[k].st_sizeZ,2*n+1);
                              }else if(zcells-(n+1)<semiSt){
                                    cell[k].st_sizeZ=MIN(cell[k].st_sizeZ,2*(zcells-(n+1))+1);
                              }
                              cell[k].st_sizeZ=MIN(cell[k].st_sizeZ,2*cell[k].distsolz-1);
                        }
				
				if(mesh->xcells<sim->order){cell[k].st_sizeX=1;}
				if(mesh->ycells<sim->order){cell[k].st_sizeY=1;}
				if(mesh->zcells<sim->order){cell[k].st_sizeZ=1;}
					
                  }
            }
	}
	//Initialization of cell stencils
	for(k=0;k<mesh->ncells;k++){
		for(p=0;p<9;p++){
			cell[k].stX[p]=-1;
			cell[k].stY[p]=-1;
                  cell[k].stZ[p]=-1;
		}
	}

      for(n=0;n<zcells;n++){
            for(m=0;m<ycells;m++){
                  for(l=0;l<xcells;l++){

                        k = l + m*xcells + n*xcells*ycells;
                        k2d=l + m*xcells;

                        //x setencils
                        for(p=0;p<cell[k].st_sizeX;p++){
                              if(mesh->periodicX==0){
                                    cell[k].stX[p]=l-((cell[k].st_sizeX-1)/2)+p;
                              }else{
                                    cell[k].stX[p]=l-((cell[k].st_sizeX-1)/2)+p;
                                    if(cell[k].stX[p]<0){
                                          cell[k].stX[p]+=xcells*(1);
                                    }
                                    if(cell[k].stX[p]>xcells-1){
                                          cell[k].stX[p]+=xcells*(-1);

                                    }

                              }
                              cell[k].stX[p]+=xcells*m;
                              cell[k].stX[p]+=xcells*ycells*n;
                        }


                        //y stencils
                        for(p=0;p<cell[k].st_sizeY;p++){
                              if(mesh->periodicY==0){
                                    cell[k].stY[p]=m-((cell[k].st_sizeY-1)/2)+p;
                              }else{
                                    cell[k].stY[p]=m-((cell[k].st_sizeY-1)/2)+p;
                                    if(cell[k].stY[p]<0){
                                          cell[k].stY[p]+=ycells*(1);
                                    }
                                    if(cell[k].stY[p]>ycells-1){
                                          cell[k].stY[p]+=ycells*(-1);

                                    }

                              }
                              cell[k].stY[p]=xcells*cell[k].stY[p]+l;
                              cell[k].stY[p]+=xcells*ycells*n;
                        }

                        //z stencils
                        for(p=0;p<cell[k].st_sizeZ;p++){
                              if(mesh->periodicZ==0){
                                    cell[k].stZ[p]=n-((cell[k].st_sizeZ-1)/2)+p;
                              }else{
                                    cell[k].stZ[p]=n-((cell[k].st_sizeZ-1)/2)+p;
                                    if(cell[k].stZ[p]<0){
                                          cell[k].stZ[p]+=zcells*(1);
                                    }
                                    if(cell[k].stZ[p]>zcells-1){
                                          cell[k].stZ[p]+=zcells*(-1);

                                    }

                              }
                              cell[k].stZ[p]=ycells*xcells*cell[k].stZ[p]+k2d;
                        }




                  }
            }

      }


      return 1;

}


int read_solids(t_mesh *mesh, t_solid *solids, const char *folder_path) {
    int i, j, k;
    FILE *fp;
    char fname[1024];
    char buffer[256];
    double aux1, aux2, aux3, aux4;
    t_triangle *triangle;
    t_stl *stl;

    snprintf(fname, sizeof(fname), "%s/solid_list.txt", folder_path);
    fp = fopen(fname, "r");

    if (fp == NULL) {
        printf("%s No solids are found in folder %s.\n", WAR, folder_path);
        solids->nsolid = 0;
        return 0;
    } else {
        if (fscanf(fp, "%*s %d", &solids->nsolid) != 1) {
            printf("Error when reading number of solids.\n");
            solids->nsolid = 0;
            fclose(fp);
            return 0;
        } else {
            printf("The number of solids domains is: %d\n", solids->nsolid);

            for (i = 0; i < solids->nsolid; i++) {
                if (fscanf(fp, "%s", buffer) != 1) {
                    printf("%s Error when reading solid names.\n", WAR);
                }
                solids->filename[i] = strdup(buffer);
                printf("Solid %d is located at %s\n", i, solids->filename[i]);
            }
        }
        fclose(fp);

        solids->stl = (t_stl *)malloc(solids->nsolid * sizeof(t_stl));
        stl = solids->stl;

        for (i = 0; i < solids->nsolid; i++) {
            snprintf(stl[i].name, sizeof(stl[i].name), "%s", solids->filename[i]);
            fp = fopen(stl[i].name, "r");
            if (fp == NULL) {
                printf("Error opening file %s.\n", solids->filename[i]);
                continue;
            }

            // Count the number of triangles
            stl[i].ntri = 0;
            while (fgets(buffer, sizeof(buffer), fp)) {
                if (strncmp(buffer, "facet normal", 12) == 0) {
                    stl[i].ntri++;
                }
            }
            
            printf("The number of triangles of solid %d is  %d\n", i, stl[i].ntri);
            
            fseek(fp, 0, SEEK_SET); // Reset file pointer to start

            stl[i].triangle = (t_triangle *)malloc(stl[i].ntri * sizeof(t_triangle));
            triangle = stl[i].triangle;

            // Initialize bounding box
            for (k = 0; k < 3; k++) {
                stl[i].Xmin[k] = 9999999999999.0;
                stl[i].Xmax[k] = -9999999999999.0;
            }

            int current_triangle = 0;
            while (fgets(buffer, sizeof(buffer), fp)) {
                if (strncmp(buffer, "facet normal", 12) == 0) {
                    // Read the normal vector
                    sscanf(buffer, "facet normal %lf %lf %lf",
                           &triangle[current_triangle].nr[0],
                           &triangle[current_triangle].nr[1],
                           &triangle[current_triangle].nr[2]);

                    // Skip the "outer loop" line
                    if (fgets(buffer, sizeof(buffer), fp) == NULL) {
                        fprintf(stderr, "Error or premature end of file when reading 'outer loop'\n");
                        break;
                    }

                    // Read the first vertex
                    if (fgets(buffer, sizeof(buffer), fp) != NULL) {
                        sscanf(buffer, " %*s %lf %lf %lf",
                               &triangle[current_triangle].p1[0],
                               &triangle[current_triangle].p1[1],
                               &triangle[current_triangle].p1[2]);
                    } else {
                        fprintf(stderr, "Error or premature end of file when reading first vertex\n");
                        break;
                    }

                    // Read the second vertex
                    if (fgets(buffer, sizeof(buffer), fp) != NULL) {
                        sscanf(buffer, " %*s %lf %lf %lf",
                               &triangle[current_triangle].p2[0],
                               &triangle[current_triangle].p2[1],
                               &triangle[current_triangle].p2[2]);
                    } else {
                        fprintf(stderr, "Error or premature end of file when reading second vertex\n");
                        break;
                    }

                    // Read the third vertex
                    if (fgets(buffer, sizeof(buffer), fp) != NULL) {
                        sscanf(buffer, " %*s %lf %lf %lf",
                               &triangle[current_triangle].p3[0],
                               &triangle[current_triangle].p3[1],
                               &triangle[current_triangle].p3[2]);
                    } else {
                        fprintf(stderr, "Error or premature end of file when reading third vertex\n");
                        break;
                    }

                    // Skip the "endloop" line
                    if (fgets(buffer, sizeof(buffer), fp) == NULL) {
                        fprintf(stderr, "Error or premature end of file when reading 'endloop'\n");
                        break;
                    }

                    // Skip the "endfacet" line
                    if (fgets(buffer, sizeof(buffer), fp) == NULL) {
                        fprintf(stderr, "Error or premature end of file when reading 'endfacet'\n");
                        break;
                    }

                    // Move to the next triangle
                    current_triangle++;
                }
            }
  

            fclose(fp);
            
            for(j=0;j<solids->stl[i].ntri;j++){
                printf("Normal %lf %lf %lf \n",triangle[j].nr[0], triangle[j].nr[1], triangle[j].nr[2]);
                printf("vertex %lf %lf %lf \n",triangle[j].p1[0], triangle[j].p1[1], triangle[j].p1[2]);
                printf("vertex %lf %lf %lf \n",triangle[j].p2[0], triangle[j].p2[1], triangle[j].p2[2]);
                printf("vertex %lf %lf %lf \n",triangle[j].p3[0], triangle[j].p3[1], triangle[j].p3[2]);
            }

            // Process triangles
            for (j = 0; j < stl[i].ntri; j++) {
                triangle[j].absnr = sqrt(triangle[j].nr[0] * triangle[j].nr[0] + triangle[j].nr[1] * triangle[j].nr[1] + triangle[j].nr[2] * triangle[j].nr[2]);
                triangle[j].outside = 0;

                  for(k=0;k<3;k++){
                        aux1=MIN(triangle[j].p1[k],triangle[j].p2[k]);
                        aux2=MIN(aux1,triangle[j].p3[k]);
                        stl[i].Xmin[k]=MIN(aux2,stl[i].Xmin[k]);

                        aux3=MAX(triangle[j].p1[k],triangle[j].p2[k]);
                        aux4=MAX(aux3,triangle[j].p3[k]);
                        stl[i].Xmax[k]=MAX(aux4,stl[i].Xmax[k]);



                        if(k==0){
                              triangle[j].imin[k]= aux2/mesh->dx;// DEBERIAMOS CONSIDERAR UN dx/2 
                              triangle[j].imax[k]= aux4/mesh->dx;
                              if(triangle[j].imax[k]-triangle[j].imin[k] < MAX(mesh->sim->order-1,1)){
                                    triangle[j].imin[k] = (aux2+aux4)/(2.0*mesh->dx) - (mesh->sim->order-1)/2;
                                    triangle[j].imax[k] = triangle[j].imin[k] + mesh->sim->order;
                                    //printf("triangle %d: dif %d   \n\n", j,triangle[j].imax[k]-triangle[j].imin[k]);
                                    //printf("triangle %d: %d %d  \n\n", j,triangle[j].imin[k],triangle[j].imax[k]);
                              }
                              triangle[j].imin[k]=MAX(triangle[j].imin[k],0);                   // to do: poner flag si el triangulo esta fuera, para sacarlo del calculo, y quitar esos min/max. Poner una minima separacion entre imin/imax
                              triangle[j].imin[k]=MIN(triangle[j].imin[k],mesh->xcells-1);
                              triangle[j].imax[k]=MAX(triangle[j].imax[k],0);
                              triangle[j].imax[k]=MIN(triangle[j].imax[k],mesh->xcells-1);
                              if(aux2<0.0||aux4>mesh->Lx){
                                    triangle[j].outside=1;
                              }
                              //printf("triangle %d: %lf %lf  \n", j,aux2,aux4);
                              //printf("triangle %d: %d %d  \n\n", j,triangle[j].imin[k],triangle[j].imax[k]);
                        }else if(k==1){
                              triangle[j].imin[k]= aux2/mesh->dy;
                              triangle[j].imax[k]= aux4/mesh->dy;
                              if(triangle[j].imax[k]-triangle[j].imin[k] < MAX(mesh->sim->order-1,1)){
                                    triangle[j].imin[k] = (aux2+aux4)/(2.0*mesh->dx) - (mesh->sim->order-1)/2;
                                    triangle[j].imax[k] = triangle[j].imin[k] + mesh->sim->order;
                              }
                              triangle[j].imin[k]=MAX(triangle[j].imin[k],0);
                              triangle[j].imin[k]=MIN(triangle[j].imin[k],mesh->ycells-1);
                              triangle[j].imax[k]=MAX(triangle[j].imax[k],0);
                              triangle[j].imax[k]=MIN(triangle[j].imax[k],mesh->ycells-1);
                              if(aux2<0.0||aux4>mesh->Ly){
                                    triangle[j].outside=1;
                              }
                              //printf("triangle %d: %lf %lf  \n", j,aux2,aux4);
                              //printf("triangle %d: %d %d  \n", j,triangle[j].imin[k],triangle[j].imax[k]);
                        }else{
                              triangle[j].imin[k]= aux2/mesh->dz;
                              triangle[j].imax[k]= aux4/mesh->dz;
                              if(triangle[j].imax[k]-triangle[j].imin[k] < MAX(mesh->sim->order-1,1)){
                                    triangle[j].imin[k] = (aux2+aux4)/(2.0*mesh->dx) - (mesh->sim->order-1)/2;
                                    triangle[j].imax[k] = triangle[j].imin[k] + mesh->sim->order;
                              }
                              triangle[j].imin[k]=MAX(triangle[j].imin[k],0);
                              triangle[j].imin[k]=MIN(triangle[j].imin[k],mesh->zcells-1);
                              triangle[j].imax[k]=MAX(triangle[j].imax[k],0);
                              triangle[j].imax[k]=MIN(triangle[j].imax[k],mesh->zcells-1);
                              if(aux2<0.0||aux4>mesh->Lz){
                                    triangle[j].outside=1;
                              }
                              //printf("triangle %d: %lf %lf  \n", j,aux2,aux4);
                              //printf("triangle %d: %d %d  \n", j,triangle[j].imin[k],triangle[j].imax[k]);
                        }
                        //triangle[j].imin[k]=(int)

                  }
                  //printf("triangle %d: %d \n", j,triangle[j].outside);

            }
            //getchar();

            for(k=0;k<3;k++){
                  if(k==0){
                        stl[i].imin[k]= stl[i].Xmin[k]/mesh->dx;
                        stl[i].imax[k]= stl[i].Xmax[k]/mesh->dx;
                        stl[i].imin[k]=MAX(stl[i].imin[k],0);
                        stl[i].imin[k]=MIN(stl[i].imin[k],mesh->xcells-1);
                        stl[i].imax[k]=MAX(stl[i].imax[k],0);
                        stl[i].imax[k]=MIN(stl[i].imax[k],mesh->xcells-1);
                        //printf("stl %d: %lf %lf  \n", i,aux2,aux4);
                        //printf("stl %d: %d %d  \n", i,stl[i].imin[k],stl[i].imax[k]);
                  }else if(k==1){
                        stl[i].imin[k]= stl[i].Xmin[k]/mesh->dy;
                        stl[i].imax[k]= stl[i].Xmax[k]/mesh->dy;
                        stl[i].imin[k]=MAX(stl[i].imin[k],0);
                        stl[i].imin[k]=MIN(stl[i].imin[k],mesh->ycells-1);
                        stl[i].imax[k]=MAX(stl[i].imax[k],0);
                        stl[i].imax[k]=MIN(stl[i].imax[k],mesh->ycells-1);
                        //printf("stl %d: %lf %lf  \n", i,aux2,aux4);
                        //printf("stl %d: %d %d  \n", i,stl[i].imin[k],stl[i].imax[k]);
                  }else{
                        stl[i].imin[k]= stl[i].Xmin[k]/mesh->dz;
                        stl[i].imax[k]= stl[i].Xmax[k]/mesh->dz;
                        stl[i].imin[k]=MAX(stl[i].imin[k],0);
                        stl[i].imin[k]=MIN(stl[i].imin[k],mesh->zcells-1);
                        stl[i].imax[k]=MAX(stl[i].imax[k],0);
                        stl[i].imax[k]=MIN(stl[i].imax[k],mesh->zcells-1);
                        //printf("stl %d: %lf %lf  \n", i,aux2,aux4);
                        //printf("stl %d: %d %d  \n", i,stl[i].imin[k],stl[i].imax[k]);
                  }
            }
            printf(" The bounding box of solid %d is: \n (x,y,z)_min=(%lf,%lf,%lf)\n (x,y,z)_max=(%lf,%lf,%lf) \n", i,solids->stl[i].Xmin[0],solids->stl[i].Xmin[1],solids->stl[i].Xmin[2],solids->stl[i].Xmax[0],solids->stl[i].Xmax[1],solids->stl[i].Xmax[2]);
            printf(" and the respective indices are:  \n (i,j,k)_min=(%d,%d,%d)\n (i,j,k)_max=(%d,%d,%d) \n", solids->stl[i].imin[0],solids->stl[i].imin[1],solids->stl[i].imin[2],solids->stl[i].imax[0],solids->stl[i].imax[1],solids->stl[i].imax[2]);
		printf("\n");
	}

      printf("%s Solid domains have been successfully read\n",OK);
    }
    return 1;
}


int read_solids_txt(t_mesh *mesh,t_solid *solids, const char *folder_path){

	int i,j,k;
	FILE *fp;
	char fname[1024];
	char buffer[256];
      double aux1,aux2,aux3,aux4;
      t_triangle *triangle;
      t_stl *stl;

	snprintf(fname, sizeof(fname), "%s/solid_list.txt", folder_path);
	fp = fopen(fname,"r");

	if (fp == NULL)
	{
		printf("%s No solids are found in folder case/.. \n",WAR);
		solids->nsolid=0;
	}else{

	if (fscanf(fp, "%*s %d", &solids->nsolid) != 1) {
		printf("Error when reading number of solids.\n");
		solids->nsolid = 0;
	}else{
		printf(" The number of solids domains is: %d \n",solids->nsolid);

		for(i=0;i<solids->nsolid;i++)
		{
		  if (fscanf(fp,"%s", buffer) != 1){printf("%s Error when reading solid names \n",WAR);}
		  solids->filename[i] = strdup(buffer);
		  printf(" Solid %d is located at %s \n",i, solids->filename[i]);
		}

	}
	printf("\n");
	fclose(fp);


	solids->stl=(t_stl*)malloc(solids->nsolid*sizeof(t_stl));
	stl=solids->stl;

	for(i=0;i<solids->nsolid;i++){
		//old: strncpy(solids->stl[i].name, solids->filename[i],256);
		snprintf(solids->stl[i].name, 256, "%s", solids->filename[i]);
		fp = fopen(solids->stl[i].name,"r");
		if (fscanf(fp, "%*s %d", &solids->stl[i].ntri) != 1) {
		printf("Error when reading triangle number %s.\n", solids->filename[i]);
		solids->stl[i].ntri = 0;
		}
		if (fscanf(fp, "%*s %d", &solids->stl[i].nver) != 1) {
		printf("Error when reading vertex number %s.\n", solids->filename[i]);
		solids->stl[i].nver = 0;
		}
		printf(" Solid %d (%s) has %d facets and %d vertices \n",i, solids->filename[i],solids->stl[i].ntri,solids->stl[i].nver);

            solids->stl[i].triangle=(t_triangle*)malloc(solids->stl[i].ntri*sizeof(t_triangle));
            triangle=solids->stl[i].triangle;


            for(k=0;k<3;k++){
                  stl[i].Xmin[k]= 9999999999999.0;
                  stl[i].Xmax[k]=-9999999999999.0;
            }


            for(j=0;j<solids->stl[i].ntri;j++){
			if (fscanf(fp, "%lf %lf %lf", &triangle[j].nr[0], &triangle[j].nr[1], &triangle[j].nr[2]) != 3) {
			    printf("Error reading the normal of triangle %d for %s.\n", j, solids->filename[i]);
			    continue; // Skip this triangle if it cannot be read completely
			}
			if (fscanf(fp, "%lf %lf %lf", &triangle[j].p1[0], &triangle[j].p1[1], &triangle[j].p1[2]) != 3) {
			    printf("Error reading node P1 of triangle %d for %s.\n", j, solids->filename[i]);
			    continue; // Skip this triangle if it cannot be read completely
			}
			if (fscanf(fp, "%lf %lf %lf", &triangle[j].p2[0], &triangle[j].p2[1], &triangle[j].p2[2]) != 3) {
			    printf("Error reading node P2 of triangle %d for %s.\n", j, solids->filename[i]);
			    continue; // Skip this triangle if it cannot be read completely
			}
			if (fscanf(fp, "%lf %lf %lf", &triangle[j].p3[0], &triangle[j].p3[1], &triangle[j].p3[2]) != 3) {
			    printf("Error reading node P3 of triangle %d for %s.\n", j, solids->filename[i]);
			    continue; // Skip this triangle if it cannot be read completely
			}
                  //printf("%lf %lf %lf\n", triangle[j].nr[0],triangle[j].nr[1],triangle[j].nr[2]);
                  //printf("%lf %lf %lf\n", triangle[j].p1[0],triangle[j].p1[1],triangle[j].p1[2]);
                  //printf("%lf %lf %lf\n", triangle[j].p2[0],triangle[j].p2[1],triangle[j].p2[2]);
                  //printf("%lf %lf %lf\n", triangle[j].p3[0],triangle[j].p3[1],triangle[j].p3[2]);
                  triangle[j].absnr=pow(triangle[j].nr[0]*triangle[j].nr[0] + triangle[j].nr[1]*triangle[j].nr[1] + triangle[j].nr[2]*triangle[j].nr[2]  ,  0.5);
                  triangle[j].outside=0;
                  for(k=0;k<3;k++){
                        aux1=MIN(triangle[j].p1[k],triangle[j].p2[k]);
                        aux2=MIN(aux1,triangle[j].p3[k]);
                        stl[i].Xmin[k]=MIN(aux2,stl[i].Xmin[k]);

                        aux3=MAX(triangle[j].p1[k],triangle[j].p2[k]);
                        aux4=MAX(aux3,triangle[j].p3[k]);
                        stl[i].Xmax[k]=MAX(aux4,stl[i].Xmax[k]);



                        if(k==0){
                              triangle[j].imin[k]= aux2/mesh->dx;// DEBERIAMOS CONSIDERAR UN dx/2 
                              triangle[j].imax[k]= aux4/mesh->dx;
                              if(triangle[j].imax[k]-triangle[j].imin[k] < MAX(mesh->sim->order-1,1)){
                                    triangle[j].imin[k] = (aux2+aux4)/(2.0*mesh->dx) - (mesh->sim->order-1)/2;
                                    triangle[j].imax[k] = triangle[j].imin[k] + mesh->sim->order;
                                    //printf("triangle %d: dif %d   \n\n", j,triangle[j].imax[k]-triangle[j].imin[k]);
                                    //printf("triangle %d: %d %d  \n\n", j,triangle[j].imin[k],triangle[j].imax[k]);
                              }
                              triangle[j].imin[k]=MAX(triangle[j].imin[k],0);                   // to do: poner flag si el triangulo esta fuera, para sacarlo del calculo, y quitar esos min/max. Poner una minima separacion entre imin/imax
                              triangle[j].imin[k]=MIN(triangle[j].imin[k],mesh->xcells-1);
                              triangle[j].imax[k]=MAX(triangle[j].imax[k],0);
                              triangle[j].imax[k]=MIN(triangle[j].imax[k],mesh->xcells-1);
                              if(aux2<0.0||aux4>mesh->Lx){
                                    triangle[j].outside=1;
                              }
                              //printf("triangle %d: %lf %lf  \n", j,aux2,aux4);
                              //printf("triangle %d: %d %d  \n\n", j,triangle[j].imin[k],triangle[j].imax[k]);
                        }else if(k==1){
                              triangle[j].imin[k]= aux2/mesh->dy;
                              triangle[j].imax[k]= aux4/mesh->dy;
                              if(triangle[j].imax[k]-triangle[j].imin[k] < MAX(mesh->sim->order-1,1)){
                                    triangle[j].imin[k] = (aux2+aux4)/(2.0*mesh->dx) - (mesh->sim->order-1)/2;
                                    triangle[j].imax[k] = triangle[j].imin[k] + mesh->sim->order;
                              }
                              triangle[j].imin[k]=MAX(triangle[j].imin[k],0);
                              triangle[j].imin[k]=MIN(triangle[j].imin[k],mesh->ycells-1);
                              triangle[j].imax[k]=MAX(triangle[j].imax[k],0);
                              triangle[j].imax[k]=MIN(triangle[j].imax[k],mesh->ycells-1);
                              if(aux2<0.0||aux4>mesh->Ly){
                                    triangle[j].outside=1;
                              }
                              //printf("triangle %d: %lf %lf  \n", j,aux2,aux4);
                              //printf("triangle %d: %d %d  \n", j,triangle[j].imin[k],triangle[j].imax[k]);
                        }else{
                              triangle[j].imin[k]= aux2/mesh->dz;
                              triangle[j].imax[k]= aux4/mesh->dz;
                              if(triangle[j].imax[k]-triangle[j].imin[k] < MAX(mesh->sim->order-1,1)){
                                    triangle[j].imin[k] = (aux2+aux4)/(2.0*mesh->dx) - (mesh->sim->order-1)/2;
                                    triangle[j].imax[k] = triangle[j].imin[k] + mesh->sim->order;
                              }
                              triangle[j].imin[k]=MAX(triangle[j].imin[k],0);
                              triangle[j].imin[k]=MIN(triangle[j].imin[k],mesh->zcells-1);
                              triangle[j].imax[k]=MAX(triangle[j].imax[k],0);
                              triangle[j].imax[k]=MIN(triangle[j].imax[k],mesh->zcells-1);
                              if(aux2<0.0||aux4>mesh->Lz){
                                    triangle[j].outside=1;
                              }
                              //printf("triangle %d: %lf %lf  \n", j,aux2,aux4);
                              //printf("triangle %d: %d %d  \n", j,triangle[j].imin[k],triangle[j].imax[k]);
                        }
                        //triangle[j].imin[k]=(int)

                  }
                  //printf("triangle %d: %d \n", j,triangle[j].outside);

            }
            //getchar();

            for(k=0;k<3;k++){
                  if(k==0){
                        stl[i].imin[k]= stl[i].Xmin[k]/mesh->dx;
                        stl[i].imax[k]= stl[i].Xmax[k]/mesh->dx;
                        stl[i].imin[k]=MAX(stl[i].imin[k],0);
                        stl[i].imin[k]=MIN(stl[i].imin[k],mesh->xcells-1);
                        stl[i].imax[k]=MAX(stl[i].imax[k],0);
                        stl[i].imax[k]=MIN(stl[i].imax[k],mesh->xcells-1);
                        //printf("stl %d: %lf %lf  \n", i,aux2,aux4);
                        //printf("stl %d: %d %d  \n", i,stl[i].imin[k],stl[i].imax[k]);
                  }else if(k==1){
                        stl[i].imin[k]= stl[i].Xmin[k]/mesh->dy;
                        stl[i].imax[k]= stl[i].Xmax[k]/mesh->dy;
                        stl[i].imin[k]=MAX(stl[i].imin[k],0);
                        stl[i].imin[k]=MIN(stl[i].imin[k],mesh->ycells-1);
                        stl[i].imax[k]=MAX(stl[i].imax[k],0);
                        stl[i].imax[k]=MIN(stl[i].imax[k],mesh->ycells-1);
                        //printf("stl %d: %lf %lf  \n", i,aux2,aux4);
                        //printf("stl %d: %d %d  \n", i,stl[i].imin[k],stl[i].imax[k]);
                  }else{
                        stl[i].imin[k]= stl[i].Xmin[k]/mesh->dz;
                        stl[i].imax[k]= stl[i].Xmax[k]/mesh->dz;
                        stl[i].imin[k]=MAX(stl[i].imin[k],0);
                        stl[i].imin[k]=MIN(stl[i].imin[k],mesh->zcells-1);
                        stl[i].imax[k]=MAX(stl[i].imax[k],0);
                        stl[i].imax[k]=MIN(stl[i].imax[k],mesh->zcells-1);
                        //printf("stl %d: %lf %lf  \n", i,aux2,aux4);
                        //printf("stl %d: %d %d  \n", i,stl[i].imin[k],stl[i].imax[k]);
                  }
            }
            printf(" The bounding box of solid %d is: \n (x,y,z)_min=(%lf,%lf,%lf)\n (x,y,z)_max=(%lf,%lf,%lf) \n", i,solids->stl[i].Xmin[0],solids->stl[i].Xmin[1],solids->stl[i].Xmin[2],solids->stl[i].Xmax[0],solids->stl[i].Xmax[1],solids->stl[i].Xmax[2]);
            printf(" and the respective indices are:  \n (i,j,k)_min=(%d,%d,%d)\n (i,j,k)_max=(%d,%d,%d) \n", solids->stl[i].imin[0],solids->stl[i].imin[1],solids->stl[i].imin[2],solids->stl[i].imax[0],solids->stl[i].imax[1],solids->stl[i].imax[2]);
		fclose(fp);
		printf("\n");
	}

      printf("%s Solid domains have been successfully read\n",OK);


	}
	return 1;
}

void set_velocity_field(t_mesh *mesh, t_sim *sim){
	/**This is only for linear transport and
	 * it only works with rho = 1.0
	 * */

	int i;
	t_wall *wall;
	t_cell *cell;
	double lambda_max,dl;
	double rho;

	rho=1.0;
	cell=mesh->cell;
	for(i=0;i<mesh->ncells;i++){
		cell->U[0]=rho;
		cell++;
	}
	lambda_max=0.0;
	wall=mesh->wall;
	for(i=0;i<mesh->nwalls;i++){
            if(wall->nx<TOL4 && wall->nz<TOL4){
            //y-wall
                  wall->fR_star[0]=(wall->cellR->U[2]+wall->cellL->U[2])*0.5;
                  wall->fL_star[0]= wall->fR_star[0];
            }else if (wall->nz<TOL4) {
                  //x-wall
                  wall->fR_star[0]=(wall->cellR->U[1]+wall->cellL->U[1])*0.5;
                  wall->fL_star[0]= wall->fR_star[0];
            }else{
                  //z-wall
                  wall->fR_star[0]=(wall->cellR->U[3]+wall->cellL->U[3])*0.5;
                  wall->fL_star[0]= wall->fR_star[0];
            }

		lambda_max=MAX(lambda_max,wall->fR_star[0]);
		wall++;
	}

	dl=MIN(mesh->dx,mesh->dy);

	sim->dt=sim->CFL*dl/lambda_max;

}

void set_velocity(t_mesh *mesh, t_sim *sim){

	int i;
	t_wall *wall;

	wall=mesh->wall;
	for(i=0;i<mesh->nwalls;i++){
            if(wall->nx<TOL4 && wall->nz<TOL4){
            //y-wall
                  wall->vel=mesh->u_y;
            }else if (wall->nz<TOL4) {
                  //x-wall
                  wall->vel=mesh->u_x;
            }else{
                  //z-wall
                  wall->vel=mesh->u_z;
            }

		wall++;
	}

}


void read_config(t_mesh *mesh, t_sim *sim, const char *folder_path){
	
	FILE *file_input;
	char fname[1024],errormsg[1024];
	
	sprintf(errormsg,"Read error in configure.input \n");
	
	snprintf(fname, sizeof(fname), "%s/configure.input", folder_path);
      file_input = fopen(fname,"r"); // "r"= only read
	if (!file_input) {
		printf("Error opening configure.input %s\n", fname);
		getchar();
		exit(1);
	}
	if (fscanf(file_input, "%*s") != 0) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &sim->tf) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &sim->tVolc) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &sim->CFL) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &sim->order) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s") != 0) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->xcells) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->ycells) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->zcells) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &mesh->Lx) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &mesh->Ly) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &mesh->Lz) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s") != 0) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->bc[0]) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->bc[1]) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->bc[2]) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->bc[3]) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->bc[4]) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %d", &mesh->bc[5]) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s") != 0) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &mesh->u_x) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &mesh->u_y) != 1) { printf("%s",errormsg); }
	if (fscanf(file_input, "%*s %lf", &mesh->u_z) != 1) { printf("%s",errormsg); }
      fclose(file_input);

}


void print_info(t_mesh *mesh, t_sim *sim, const char *folder_path){

      printf("\n\e[94m Authors:\n  - Adrián Navas Montilla\n  - Isabel Echeverribar \n");
      printf(" Copyright (C) 2019-2024 The authors.   \n\n");
      //printf("License type: The 3-Clause BSD License \nRedistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: \n 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. \n 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. \n 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.\n\n");
      //printf("This software is provided by the copyright holders and contributors “as is” and any express or implied warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose are disclaimed. In no event shall the copyright holder or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.\e[0m\n");

      printf(" \n");
      printf(" \e[4mSIMULATION SETUP:\e[0m\n");
	printf(" Folder path: %s \n",folder_path);
#if TYPE_REC==0
	printf(" WENO reconstruction of order %d is chosen. \n",sim->order);
#endif
#if TYPE_REC==1
	printf(" TENO reconstruction of order %d is chosen. \n",sim->order);
#endif
#if TYPE_REC==2
	printf("%s UWC (optimal weights) reconstruction of order %d is chosen. \n",WAR,sim->order);
#endif
      printf(" Final time: %lf\n",sim->tf);
      printf(" CFL: %lf\n",sim->CFL);
      printf(" Number of cells X: %d\n",mesh->xcells);
      printf(" Number of cells Y: %d\n",mesh->ycells);
      printf(" Number of cells Z: %d\n",mesh->zcells);
      printf(" Domain size: %lf x %lf x %lf \n",mesh->Lx,mesh->Ly,mesh->Lz);
      printf(" Boundaries (1: periodic, 2: user defined, 3: transmissive, 4: solid wall): \n"); // 1: periodic, 2: user defined, 3: transmissive, 4: solid wall
      printf(" Face_1(-y): %d \n",mesh->bc[0]);
      printf(" Face_2(+x): %d \n",mesh->bc[1]);
      printf(" Face_3(+y): %d \n",mesh->bc[2]);
      printf(" Face_4(-x): %d \n",mesh->bc[3]);
      printf(" Face_5(-z): %d \n",mesh->bc[4]);
      printf(" Face_6(+z): %d \n",mesh->bc[5]);
#if EQUATION_SYSTEM == 0
      printf("%s LINEAR TRANPORT IS ACTIVE. \n",WAR);
      printf(" Linear transport velocity: \n");
      printf(" u_x: %lf \n",mesh->u_x);
      printf(" u_y: %lf \n",mesh->u_y);
      printf(" u_z: %lf \n",mesh->u_z);
#endif

      printf(" \n");

      printf("%s Configuration file has been read \n",OK);

#if ST!=1&&SOLVER==2
	printf("%s HLLS solver cannot be used when ST=0, ST=2 or ST=3. Please use HLL or HLLC. Press any key to exit... \n",ERR);
	getchar();
	exit(1);
#endif

#if SOLVER>2
	printf("%s The solver is not selected adequately. Press any key to exit... \n",ERR);
	getchar();
	exit(1);
#endif

#if MULTICOMPONENT>0&&SOLVER==1
	printf("%s HLLC solver cannot handle multicomponent flow. Please use HLL (SOLVER = 0). Press any key to exit... \n",ERR);
	getchar();
	exit(1);
#endif

      if((mesh->bc[1]==1 && mesh->bc[3]!=1)||(mesh->bc[1]!=1 && mesh->bc[3]==1)){
            printf("%s Cyclic BC in X not properly set, only one of the boundaries is set as cyclic. The program will close when pressing a key. \n",ERR);
            getchar();
            exit(1);
	}
      if((mesh->bc[0]==1 && mesh->bc[2]!=1)||(mesh->bc[0]!=1 && mesh->bc[2]==1)){
            printf("%s Cyclic BC in Y not properly set, only one of the boundaries is set as cyclic. The program will close when pressing a key. \n",ERR);
            getchar();
            exit(1);
	}
      if((mesh->bc[4]==1 && mesh->bc[5]!=1)||(mesh->bc[4]!=1 && mesh->bc[5]==1)){
            printf("%s Cyclic BC in Z not properly set, only one of the boundaries is set as cyclic. The program will close when pressing a key. \n",ERR);
            getchar();
            exit(1);
	}


      if(mesh->xcells<=(sim->order-1)/2){
            printf("%s The number of cells in X is lower than half-the order of accuracy (xcells < (order-1)/2). If you want a 1D or 2D simulation, please be sure that xcells=1. \n",WAR);
		if(mesh->bc[1]==1){
		printf("%s The number of cells in X is too small for periodic BC. Transmissive BC are considered instead.  \n",WAR);
		mesh->bc[1]=3;
            mesh->bc[3]=3;}

	}
      if(mesh->ycells<=(sim->order-1)/2){
            printf("%s The number of cells in Y is lower than half-the order of accuracy (ycells < (order-1)/2). If you want a 1D or 2D simulation, please be sure that ycells=1. \n",WAR);
		if(mesh->bc[0]==1){
		printf("%s The number of cells in Y is too small for periodic BC. Transmissive BC are considered instead.  \n",WAR);
		mesh->bc[0]=3;
            mesh->bc[2]=3;}

	}
      if(mesh->zcells<=(sim->order-1)/2){
            printf("%s The number of cells in Z is lower than half-the order of accuracy (zcells <  (order-1)/2). If you want a 1D or 2D simulation, please be sure that zcells=1. \n",WAR);
		if(mesh->bc[4]==1){
		printf("%s The number of cells in Z is too small for periodic BC. Transmissive BC are considered instead.  \n",WAR);
		mesh->bc[4]=3;
            mesh->bc[5]=3;}

	}
    
    if((sim->order==2)||(sim->order==4)||(sim->order==6)){
            printf("%s Only odd orders are allowed (Order = 1, 3, 5 or 7). The program will close when pressing a key. \n",ERR);
            getchar();
            exit(1);
	}
	
	

}

