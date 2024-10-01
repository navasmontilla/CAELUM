/*
Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - ibmutils.c

Content:
  -This file contains all the functions related with the immersed boundary method
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
#include "ibmutils.h"



int assign_image_cells(t_mesh *mesh,t_solid *solids){
	t_cell *cell;
	int n,q;
	int imin,imax,jmin,jmax,kmin,kmax;
      double di[8],li[8];
      double aux1,aux2,aux3,sum;

	if(solids->nsolid<1){
            printf("%s In function assign_image_cells() no solids are considered\n",WAR);
      }else{

		cell=mesh->cell;
		for(n=0;n<mesh->ncells;n++){
			  if (cell[n].ghost==1) {

				if (cell[n].xim>0.0 && cell[n].xim<mesh->Lx && cell[n].yim>0.0 && cell[n].yim<mesh->Ly && cell[n].zim>0.0 && cell[n].zim<mesh->Lz) {  // We first check that image points are inside the computational domain. Habria que usar +- dx para afinar

					// The indices of the 8 closest cells to an image point are defined
					imin= MAX((cell[n].xim-mesh->dx/2.0)/mesh->dx,0);
					imax= MIN(imin+1,mesh->xcells-1);
					jmin= MAX((cell[n].yim-mesh->dy/2.0)/mesh->dy,0);
					jmax= MIN(jmin+1,mesh->ycells-1);
					kmin= MAX((cell[n].zim-mesh->dz/2.0)/mesh->dz,0);
					kmax= MIN(kmin+1,mesh->zcells-1);



					// The cell numbers of those 8 neighbors are defined
					cell[n].ni[0]= imin + jmin*mesh->xcells  + kmin*mesh->xcells*mesh->ycells;
					cell[n].ni[1]= imax + jmin*mesh->xcells  + kmin*mesh->xcells*mesh->ycells;
					cell[n].ni[2]= imax + jmax*mesh->xcells  + kmin*mesh->xcells*mesh->ycells;
					cell[n].ni[3]= imin + jmax*mesh->xcells  + kmin*mesh->xcells*mesh->ycells;
					cell[n].ni[4]= imin + jmin*mesh->xcells  + kmax*mesh->xcells*mesh->ycells;
					cell[n].ni[5]= imax + jmin*mesh->xcells  + kmax*mesh->xcells*mesh->ycells;
					cell[n].ni[6]= imax + jmax*mesh->xcells  + kmax*mesh->xcells*mesh->ycells;
					cell[n].ni[7]= imin + jmax*mesh->xcells  + kmax*mesh->xcells*mesh->ycells;


					// The weights for the interpolation at the image points are computed
					sum=0.0;
					for(q=0;q<8;q++){ // loop over neighbor cells
						aux1=cell[n].xim-cell[cell[n].ni[q]].xc;
						aux2=cell[n].yim-cell[cell[n].ni[q]].yc;
						aux3=cell[n].zim-cell[cell[n].ni[q]].zc;

						//if(n==594473){
						//      printf("cell image : %lf %lf %lf\n", cell[n].xim,cell[n].yim,cell[n].zim );
						//      printf("aux1 %lf aux2 %lf aux3 %lf\n", aux1,aux2,aux3);
						//      getchar();
						//}

						di[q]=pow( aux1*aux1+aux2*aux2+aux3*aux3, 0.5 ); // this is the distance between the image point and each neighbor cell
						if(cell[cell[n].ni[q]].ghost!=1){
							li[q]=1.0/(di[q]*di[q]+TOL14);
						}else{
							li[q]=0.0;
						}
						sum = sum+li[q];

					}

					if(sum<TOL14){
						cell[n].type=0;
						cell[n].ghost=0;
					}else{
						for(q=0;q<8;q++){
							cell[n].li[q]=li[q]/sum; // weighting coefficients
						}
					}




				}else{
				//if(cell[n].zc<mesh->dz){
				//      printf("cell %d,: %d %d %d \n", n,cell[n].l,cell[n].m,cell[n].n );
				//      printf("cell image : %lf %lf %lf\n", cell[n].xim,cell[n].yim,cell[n].zim );
				//	}
					cell[n].type=0;
					cell[n].ghost=0;
				}

			  }
		}

		//printf("%s Image points have been identified \n\n",OK);
	}

	return 1;
}


int update_ghost_cells(t_sim *sim,t_mesh *mesh,t_solid *solids){
	t_cell *cell;
      t_triangle *triangle;
	int n,k,q;
      double auxval[sim->nvar];
      double dotprod;

	if(solids->nsolid>0){
		cell=mesh->cell;
		for(n=0;n<mesh->ncells;n++){
			  if (cell[n].ghost==1) {
				triangle=cell[n].tri;  // triangular facet associated to a ghost cell
				for(k=0;k<sim->nvar;k++){
					auxval[k]=0.0;
					for(q=0;q<8;q++){
						auxval[k]= auxval[k] + cell[n].li[q]* cell[cell[n].ni[q]].U[k]; //interpolated variables at image point
						//printf("li:%lf U:%lf \n", cell[n].li[q],cell[cell[n].ni[q]].U[k] );
					}
					//getchar();
				}


				dotprod=triangle->nr[0]*auxval[1]+triangle->nr[1]*auxval[2]+triangle->nr[2]*auxval[3];
				for(k=0;k<sim->nvar;k++){
					if(k==1||k==2||k==3){
						cell[n].U[k]=auxval[k]-2.0*dotprod*triangle->nr[k-1]; //this is a reflection for vector variables u_r=u-2*(u·n)n, which allows to impose the Dirichlet BC of zero velocity at solid faces
					}else{
						cell[n].U[k]=auxval[k]; //non-vector variables are assigned equal
					}

				}



			  }
		}

		//printf("%s Ghost cell values have been computed \n\n",OK);

      }

	return 1;
}


int update_wall_type(t_mesh *mesh,t_solid *solids){
	t_wall *wall;
	int n;

	if(solids->nsolid<1){
            printf("%s In function update_wall_type() no solids are considered\n",WAR);
      }else{
		for(n=0;n<mesh->nwalls;n++){
			wall=&(mesh->wall[n]);

			if(wall->cellL->ghost>0 && wall->cellR->ghost>0){       //left and right ghost
				wall->wtype=0;
			}else if(wall->cellR->type==0 && wall->cellL->ghost>0){ //left ghost and right solid
				wall->wtype=0;
			}else if(wall->cellL->type==0 && wall->cellR->ghost>0){ //left solid and right ghost
				wall->wtype=0;
			}else if(wall->cellL->type==0 && wall->cellR->type==0){ //left and right solids
				wall->wtype=0;
			}


		}

	}

      return 1;

}

