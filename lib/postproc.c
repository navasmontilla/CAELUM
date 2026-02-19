/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - postproc.c

Content:
  -This file contains all the post-processing utilities
  
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
#include "postproc.h"
#include "closures.h"
#include "numcore.h"


int write_vtk(t_mesh *mesh, char *filename){


	int i,j;
	FILE *fp;
#if EQUATION_SYSTEM == 2 
	double gamma,theta,u,v,w;
#endif
	fp=fopen(filename,"w");

	// Write file header
	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"Output file from Euler\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

	// Write node coordinates
	fprintf(fp,"POINTS %d double \n",mesh->nnodes);
	for (i=0;i<mesh->nnodes;i++){
	   fprintf(fp,"%lf %lf %lf\n", mesh->node[i].x, mesh->node[i].y, mesh->node[i].z);
	}

	//Write cell-node connectivity
	fprintf(fp,"CELLS %d %d \n",mesh->ncells,mesh->ncells*(8+1));
	for(i=0;i<mesh->ncells;i++){
		fprintf(fp,"8 %d %d %d %d %d %d %d %d\n",mesh->cell[i].n1,mesh->cell[i].n2,mesh->cell[i].n3,mesh->cell[i].n4,mesh->cell[i].n5,mesh->cell[i].n6,mesh->cell[i].n7,mesh->cell[i].n8);
	}
	////////////////////////////////////////////////////
	////////////// C O N E C T I V I T Y ///////////////
	////////////////////////////////////////////////////

	// Cell types
	fprintf(fp,"CELL_TYPES %d \n",mesh->ncells);
	for(i=0;i<mesh->ncells;i++){
	   fprintf(fp,"%d \n",12);
	}

	//Cell data
	fprintf(fp,"CELL_DATA %d \n",mesh->ncells);


#if EQUATION_SYSTEM == 2 
      #if print_RHO
	fprintf(fp,"SCALARS rho DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e\n",mesh->cell[j].U[0]);
	}
      #endif

	for(j=0;j<mesh->ncells;j++){
		#if MULTICOMPONENT
			#if MULTI_TYPE==1
				gamma=mesh->cell[j].U[5]/mesh->cell[j].U[0];
			#else
				gamma=1.0+1.0/(mesh->cell[j].U[5]/mesh->cell[j].U[0]);
			#endif
		#else
			gamma=_gamma_;
		#endif
		u=mesh->cell[j].U[1]/mesh->cell[j].U[0];
		v=mesh->cell[j].U[2]/mesh->cell[j].U[0];
		w=mesh->cell[j].U[3]/mesh->cell[j].U[0];
            mesh->cell[j].pres=pressure_from_energy(gamma, mesh->cell[j].U[4], u, v, w, mesh->cell[j].U[0], mesh->cell[j].zc);
		//mesh->cell[j].pres=(gamma-1.0)*(mesh->cell[j].U[4]-0.5*mesh->cell[j].U[0]*(mesh->cell[j].U[1]*mesh->cell[j].U[1]+mesh->cell[j].U[2]*mesh->cell[j].U[2]+mesh->cell[j].U[3]*mesh->cell[j].U[3])/(mesh->cell[j].U[0]*mesh->cell[j].U[0]));
      }
      #if print_PRESSURE
	fprintf(fp,"SCALARS pres DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
      for(j=0;j<mesh->ncells;j++){
		fprintf(fp,"%14.14e \n",mesh->cell[j].pres);
      }
      #endif

      #if print_OVERPRESSURE
	fprintf(fp,"SCALARS d_pres DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e \n",mesh->cell[j].pres - mesh->cell[j].prese);
	}
      #endif

      #if print_VELOCITY
	fprintf(fp,"SCALARS U DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e \n",mesh->cell[j].U[1]/mesh->cell[j].U[0]);
	}

	fprintf(fp,"SCALARS V DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e \n",mesh->cell[j].U[2]/mesh->cell[j].U[0]);
	}

	fprintf(fp,"SCALARS W DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e \n",mesh->cell[j].U[3]/mesh->cell[j].U[0]);
	}
      #endif

      #if print_ENERGY
	fprintf(fp,"SCALARS E DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e \n",mesh->cell[j].U[4]);
	}
      #endif

      #if print_SOLUTES
	fprintf(fp,"SCALARS phi DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e \n",mesh->cell[j].U[5]);
	}
      #endif

      #if print_POTENTIALTEM
	fprintf(fp,"SCALARS theta DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
            theta=mesh->cell[j].pres/(_R_*mesh->cell[j].U[0])/( pow((mesh->cell[j].pres/_p0_),((_gamma_-1.0)/_gamma_)) );
	   	fprintf(fp,"%14.14e \n",theta);
	}
      #endif


#else
	fprintf(fp,"SCALARS U DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%14.14e\n",mesh->cell[j].U[0]);
	}
#endif

	fclose(fp);
	printf("%s A VTK has been written: %s\n",OK,filename);


	return 1;
}


int write_list(t_mesh *mesh, char *filename){

#if EQUATION_SYSTEM == 2 

	int l,m,n,k;
	double u,v,w,p,rho,phi,gamma,theta;
	FILE *fp;
	fp=fopen(filename,"w");
	// Write file header
	fprintf(fp,"VARIABLES = X, Y, Z, u, v, w, rho, p, phi, theta \n");
	fprintf(fp,"CELLS = %d, %d, %d,\n",mesh->xcells,mesh->ycells,mesh->zcells);

    for(l=0;l<mesh->xcells;l++){
		for(m=0;m<mesh->ycells;m++){
			for(n=0;n<mesh->zcells;n++){
			k = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
			u=mesh->cell[k].U[1]/mesh->cell[k].U[0];
			v=mesh->cell[k].U[2]/mesh->cell[k].U[0];
			w=mesh->cell[k].U[3]/mesh->cell[k].U[0];
			rho=mesh->cell[k].U[0];
			phi=mesh->cell[k].U[5]/mesh->cell[k].U[0];
			#if MULTICOMPONENT
				#if MULTI_TYPE==1
					gamma=phi;
				#else
					gamma=1.0+1.0/phi;
				#endif
			#else
				gamma=_gamma_;
			#endif
			p=pressure_from_energy(gamma, mesh->cell[k].U[4], u, v, w, mesh->cell[k].U[0], mesh->cell[k].zc);
			theta=p/(_R_*mesh->cell[k].U[0])/( pow((p/_p0_),((_gamma_-1.0)/_gamma_)) );
                  fprintf(fp,"%14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e\n",mesh->cell[k].xc,mesh->cell[k].yc,mesh->cell[k].zc,u,v,w,rho,p,phi,theta);
            }
		}
	}
	fclose(fp);
	
#else
	int l,m,n,k;
	FILE *fp;
	fp=fopen(filename,"w");

	// Write file header
	fprintf(fp,"VARIABLES = X, Y, Z, U \n");
	fprintf(fp,"CELLS = %d, %d, %d,\n",mesh->xcells,mesh->ycells,mesh->zcells);

    for(l=0;l<mesh->xcells;l++){
		for(m=0;m<mesh->ycells;m++){
			for(n=0;n<mesh->zcells;n++){
			k = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
			fprintf(fp,"%14.14e %14.14e %14.14e %14.14e \n",mesh->cell[k].xc,mesh->cell[k].yc,mesh->cell[k].zc,mesh->cell[k].U[0]);
            }
		}
	}
	fclose(fp);
	
#endif

	
	printf("%s A *.out file has been written: %s\n",OK,filename);

	return 1;
}

int write_list_eq(t_mesh *mesh, char *filename){
#if (EQUATION_SYSTEM == 2  && ST!=0)
	int l,m,n,k;
	double u,v,w,p,rho,phi,gamma,theta;
	FILE *fp;
	fp=fopen(filename,"w");

	// Write file header
	fprintf(fp,"VARIABLES = X, Y, Z, u, v, w, rho, p, phi, theta \n");
	fprintf(fp,"CELLS = %d, %d, %d,\n",mesh->xcells,mesh->ycells,mesh->zcells);

    for(l=0;l<mesh->xcells;l++){
		for(m=0;m<mesh->ycells;m++){
			for(n=0;n<mesh->zcells;n++){
			k = l + m*mesh->xcells + n*mesh->xcells*mesh->ycells;
			u=mesh->cell[k].Ue[1]/mesh->cell[k].Ue[0];
			v=mesh->cell[k].Ue[2]/mesh->cell[k].Ue[0];
			w=mesh->cell[k].Ue[3]/mesh->cell[k].Ue[0];
			rho=mesh->cell[k].Ue[0];
			phi=mesh->cell[k].Ue[5]/mesh->cell[k].Ue[0];
			#if MULTICOMPONENT
				#if MULTI_TYPE==1
					gamma=phi;
				#else
					gamma=1.0+1.0/phi;
				#endif
			#else
				gamma=_gamma_;
			#endif
			p=pressure_from_energy(gamma, mesh->cell[k].Ue[4], u, v, w, mesh->cell[k].Ue[0], mesh->cell[k].zc);
			theta=p/(_R_*mesh->cell[k].Ue[0])/( pow((p/_p0_),((_gamma_-1.0)/_gamma_)) );
                  fprintf(fp,"%14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e %14.14e \n",mesh->cell[k].xc,mesh->cell[k].yc,mesh->cell[k].zc,u,v,w,rho,p,phi,theta);
            }
		}
	}

	fclose(fp);
#endif
	return 1;
}




int write_geo_vtk(t_mesh *mesh, char *filename){

	int i,j;
	FILE *fp;
	fp=fopen(filename,"w");

	// Write file header
	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"Output file from Euler\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

	// Write node coordinates
	fprintf(fp,"POINTS %d double \n",mesh->nnodes);
	for (i=0;i<mesh->nnodes;i++){
	   fprintf(fp,"%lf %lf %lf\n", mesh->node[i].x, mesh->node[i].y, mesh->node[i].z);
	}

	//Write cell-node connectivity
	fprintf(fp,"CELLS %d %d \n",mesh->ncells,mesh->ncells*(8+1));
	for(i=0;i<mesh->ncells;i++){
		fprintf(fp,"8 %d %d %d %d %d %d %d %d\n",mesh->cell[i].n1,mesh->cell[i].n2,mesh->cell[i].n3,mesh->cell[i].n4,mesh->cell[i].n5,mesh->cell[i].n6,mesh->cell[i].n7,mesh->cell[i].n8);
	}
	////////////////////////////////////////////////////
	////////////// C O N E C T I V I T Y ///////////////
	////////////////////////////////////////////////////

	// Cell types
	fprintf(fp,"CELL_TYPES %d \n",mesh->ncells);
      	for(i=0;i<mesh->ncells;i++){
	   fprintf(fp,"%d \n",12);
	}

	//Cell data
	fprintf(fp,"CELL_DATA %d \n",mesh->ncells);

	//Water depth (h) difference

	fprintf(fp,"SCALARS stX INTEGER \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%d \n",mesh->cell[j].st_sizeX);
	}
	fprintf(fp,"SCALARS stY INTEGER \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%d \n",mesh->cell[j].st_sizeY);
	}
      fprintf(fp,"SCALARS stZ INTEGER \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%d \n",mesh->cell[j].st_sizeZ);
	}
      fprintf(fp,"SCALARS cellType INTEGER \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%d \n",mesh->cell[j].type);
	}
      #if ALLOW_SOLIDS==3
      fprintf(fp,"VECTORS normal double\n");
      for (j = 0; j < mesh->ncells; j++) {
          fprintf(fp, "%lf %lf %lf\n",
                  mesh->cell[j].nr[0],
                  mesh->cell[j].nr[1],
                  mesh->cell[j].nr[2]);
      }
      fprintf(fp,"SCALARS SignedDistance DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%f \n",mesh->cell[j].sdf);
	}
      #endif
      


	fclose(fp);
	printf("%s A VTK has been written: %s\n",OK,filename);


	return 1;
}


int write_extrema(t_mesh *mesh, char *filename){

#if EQUATION_SYSTEM == 2 
	int i,j;
	FILE *fp;
	fp=fopen(filename,"w");

	// Write file header
	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"Output file from Euler\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

	// Write node coordinates
	fprintf(fp,"POINTS %d double \n",mesh->nnodes);
	for (i=0;i<mesh->nnodes;i++){
	   fprintf(fp,"%lf %lf %lf\n", mesh->node[i].x, mesh->node[i].y, mesh->node[i].z);
	}

	//Write cell-node connectivity
	fprintf(fp,"CELLS %d %d \n",mesh->ncells,mesh->ncells*(8+1));
	for(i=0;i<mesh->ncells;i++){
		fprintf(fp,"8 %d %d %d %d %d %d %d %d\n",mesh->cell[i].n1,mesh->cell[i].n2,mesh->cell[i].n3,mesh->cell[i].n4,mesh->cell[i].n5,mesh->cell[i].n6,mesh->cell[i].n7,mesh->cell[i].n8);
	}
	////////////////////////////////////////////////////
	////////////// C O N E C T I V I T Y ///////////////
	////////////////////////////////////////////////////

	// Cell types
	fprintf(fp,"CELL_TYPES %d \n",mesh->ncells);
      	for(i=0;i<mesh->ncells;i++){
	   fprintf(fp,"%d \n",12);
	}

	//Cell data
	fprintf(fp,"CELL_DATA %d \n",mesh->ncells);

	fprintf(fp,"SCALARS pmax DOUBLE \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
      for(j=0;j<mesh->ncells;j++){
		fprintf(fp,"%14.14e \n",mesh->cell[j].pmax);
      }
      
      fprintf(fp,"SCALARS troubled INTEGER \n");
	fprintf(fp,"LOOKUP_TABLE DEFAULT \n");
	for(j=0;j<mesh->ncells;j++){
	   	fprintf(fp,"%d \n",mesh->cell[j].troubled);
	}

	fclose(fp);
	printf("%s A VTK has been written: %s\n",OK,filename);

#endif

	return 1;
}



void screen_info(t_mesh *mesh, t_sim *sim){
	
	printf("\n");
	printf(" T= %lf, dt= %lf\n",sim->t,sim->dt);
	mass_calculation(mesh,sim);
	energy_calculation(mesh,sim);
	printf(" Total mass: M= %14.14e\n",mesh->mass);
#if EQUATION_SYSTEM == 2 
	printf(" Total energy: E= %14.14e\n",mesh->energy);
#endif
	printf(" Mass error: (M-M0)/M0 = %14.14e\n",(mesh->mass-mesh->mass0)/mesh->mass0);
#if EQUATION_SYSTEM == 2 
	printf(" Energy error: (E-E0)/E0 = %14.14e\n",(mesh->energy-mesh->energy0)/mesh->energy0);
#endif
	printf("\n");
	
}