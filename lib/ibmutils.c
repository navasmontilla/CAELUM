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



int assign_image_cells(t_mesh *mesh){    
 
      #if (ALLOW_SOLIDS==1)||(ALLOW_SOLIDS==3)
            
            t_cell *cell;
            int n,q;
            int imin,imax,jmin,jmax,kmin,kmax;
            double di[8],li[8];
            double aux1,aux2,aux3,sum;
      
		cell=mesh->cell;
		for(n=0;n<mesh->ncells;n++){
			  if (cell[n].type==2) {

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

						di[q]=pow( aux1*aux1+aux2*aux2+aux3*aux3, 0.5 ); // this is the distance between the image point and each neighbor cell
						if(cell[cell[n].ni[q]].type==1){
							li[q]=1.0/(di[q]*di[q]+TOL14);
						}else{
							li[q]=0.0;
						}
						sum = sum+li[q];

					}

					if(sum<TOL14){
						cell[n].type=0;
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
				}

			  }
		}
            
            for(n=0;n<mesh->ncells;n++){
                  if (cell[n].type==0) {
                        cell[n].U[0]=-1.0;
                        cell[n].U[1]=0.0;
				cell[n].U[2]=0.0;
				cell[n].U[3]=0.0;
				cell[n].U[4]=0.0;
				cell[n].U[5]=0.0;
                  }
            }

      #endif
      
	return 1;
}


int update_ghost_cells(t_sim *sim,t_mesh *mesh){

      #if (ALLOW_SOLIDS==1)||(ALLOW_SOLIDS==3)
            
      	t_cell *cell;
            #if ALLOW_SOLIDS==1
            t_triangle *triangle;
            #endif
            int n,k,q;
            double auxval[sim->nvar];
            double dotprod;
      
		cell=mesh->cell;
		for(n=0;n<mesh->ncells;n++){
			  if (cell[n].type==2) {
				for(k=0;k<sim->nvar;k++){
					auxval[k]=0.0;
					for(q=0;q<8;q++){
						auxval[k]= auxval[k] + cell[n].li[q]* cell[cell[n].ni[q]].U[k]; //interpolated variables at image point
						//if (isnan(cell[n].li[q]) || isnan(cell[cell[n].ni[q]].U[k] )) {
                                    //printf("li:%lf U:%lf \n", cell[n].li[q],cell[cell[n].ni[q]].U[k] );}
					}
                              if(k==1||k==2||k==3){
                                    auxval[k]=auxval[k]/auxval[0];
                              }
					//getchar();
				}

                        #if ALLOW_SOLIDS==1
                        triangle=cell[n].tri;  // triangular facet associated to a ghost cell
				dotprod=triangle->nr[0]*auxval[1]+triangle->nr[1]*auxval[2]+triangle->nr[2]*auxval[3];
				for(k=0;k<sim->nvar;k++){
					if(k==1||k==2||k==3){
						cell[n].U[k]=auxval[k]-2.0*dotprod*triangle->nr[k-1]; //this is a reflection for vector variables u_r=u-2*(u·n)n, which allows to impose the Dirichlet BC of zero velocity at solid faces
                                    cell[n].U[k]=cell[n].U[k]*cell[n].U[0];
                              }else{
						cell[n].U[k]=auxval[k]; //non-vector variables are assigned equal
					}

				}
                        #endif
                        
                        #if ALLOW_SOLIDS==3
				dotprod=cell[n].nr[0]*auxval[1]+cell[n].nr[1]*auxval[2]+cell[n].nr[2]*auxval[3];
				for(k=0;k<sim->nvar;k++){
					if(k==1||k==2||k==3){
						cell[n].U[k]=auxval[k]-2.0*dotprod*cell[n].nr[k-1]; //this is a reflection for vector variables u_r=u-2*(u·n)n, which allows to impose the Dirichlet BC of zero velocity at solid faces
                                    cell[n].U[k]=cell[n].U[k]*cell[n].U[0];
                              }else{
						cell[n].U[k]=auxval[k]; //non-vector variables are assigned equal
					}

				}
                        #endif



			  }
		}
      
      #endif
      
	return 1;
}


int update_wall_type(t_mesh *mesh){

      #if (ALLOW_SOLIDS==1)||(ALLOW_SOLIDS==3)
            
      	t_wall *wall;
            int n;
      
		for(n=0;n<mesh->nwalls;n++){
			wall=&(mesh->wall[n]);

			if(wall->cellL->type==2 && wall->cellR->type==2){       //left and right ghost
				wall->wtype=0;
			}else if(wall->cellR->type==0 && wall->cellL->type==2){ //left ghost and right solid
				wall->wtype=0;
			}else if(wall->cellL->type==0 && wall->cellR->type==2){ //left solid and right ghost
				wall->wtype=0;
			}else if(wall->cellL->type==0 && wall->cellR->type==0){ //left and right solids
				wall->wtype=0;
			}

                  if(wall->boundary != 1){
                        if(wall->nx>TOL4){ //e imponer que no sea pared del contorno
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
                  }

		}

	#endif

      return 1;

}

