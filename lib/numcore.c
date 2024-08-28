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
  - numcore.c

Content:
  -This code file contains all the functions related with the core of the numerical scheme
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
#include "numcore.h"
#include "reconst.h"
#include "solvers.h"
#include "closures.h"
#include "ibmutils.h"


void update_cell(t_mesh *mesh, t_sim *sim){

	int i,k;
	t_cell *cell;

	cell=mesh->cell;
	for(i=0;i<mesh->ncells;i++){
            if(cell->type!=0&&cell->ghost!=1){
		for(k=0;k<sim->nvar;k++){
			cell->U[k]-=sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy + (cell->w6->fL_star[k]-cell->w5->fR_star[k])/cell->dz );
		}
            //printf("%lf\n",cell->w2->fL_star[5]);
            if (cell->U[0]<TOL14){
            	printf("celda %d: rHO: %lf \n",i,cell->U[0]);
		getchar();
            }
            }

		cell++;
	}
}


void update_cellK1(t_mesh *mesh, t_sim *sim){

	int i,k;
	t_cell *cell;

	//cell=mesh->cell;
#pragma omp parallel for default(none) private(k,cell) shared(sim,mesh)
	for(i=0;i<mesh->ncells;i++){
		cell=&(mesh->cell[i]);
            if(cell->type!=0&&cell->ghost!=1){
		for(k=0;k<sim->nvar;k++){
			cell->U_aux[k]=cell->U[k];
			cell->U[k]-=sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy + (cell->w6->fL_star[k]-cell->w5->fR_star[k])/cell->dz - cell->S[k]);
		}
            }
	}
}


void update_cellK2(t_mesh *mesh, t_sim *sim){

	int i,k;
	t_cell *cell;

	//cell=mesh->cell;
#pragma omp parallel for default(none) private(k,cell) shared(sim,mesh)
	for(i=0;i<mesh->ncells;i++){
            cell=&(mesh->cell[i]);
            if(cell->type!=0&&cell->ghost!=1){
		for(k=0;k<sim->nvar;k++){
			cell->U[k]=0.75*cell->U_aux[k]+0.25*cell->U[k]-0.25*sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy + (cell->w6->fL_star[k]-cell->w5->fR_star[k])/cell->dz - cell->S[k]);

		}
            }
	}
}


void update_cellK3(t_mesh *mesh, t_sim *sim){

	int i,k;
	t_cell *cell;

	//cell=mesh->cell;
#pragma omp parallel for default(none) private(k,cell) shared(sim,mesh)
	for(i=0;i<mesh->ncells;i++){
            cell=&(mesh->cell[i]);
            if(cell->type!=0&&cell->ghost!=1){
		for(k=0;k<sim->nvar;k++){
			cell->U[k]=(1.0/3.0)*cell->U_aux[k]+(2.0/3.0)*cell->U[k]-(2.0/3.0)*sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy + (cell->w6->fL_star[k]-cell->w5->fR_star[k])/cell->dz - cell->S[k]);

		}
		}
	}
}


int equilibrium_reconstruction(t_mesh *mesh, t_sim *sim){

	double phi3[3],phi5[5],phi7[7]; //auxiliar arrays for the weno reconstruction
	double uL,uR,vL,vR,wL,wR;
	int order;
	int n,i,j,k;
	int st[9]; //local stencil array
	t_wall *wall;
      t_cell *cell;


#pragma omp parallel for default(none) private(wall,phi3,phi5,phi7,uL,uR,vL,vR,wL,wR,order,i,j,k,st) shared(sim,mesh)
	for(n=0;n<mesh->nwalls;n++){
		wall=&(mesh->wall[n]);

            if(wall->wtype!=0){

            //RIGHT RECONSTRUCTION
            if(wall->nx<TOL4 && wall->nz<TOL4){
                  //y-wall
                  order=wall->cellR->st_sizeY;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellR->stY[j];
                  }
            }else if (wall->nz<TOL4) {
                  //x-wall
                  order=wall->cellR->st_sizeX;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellR->stX[j];
                  }
            }else{
                  //z-wall
                  order=wall->cellR->st_sizeZ;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellR->stZ[j];
                  }
            }
            if(order==1){

                  for(k=0;k<sim->nvar;k++){
                        wall->URe[k]=wall->cellR->Ue[k];
                  }

            }else if(order==3){

                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi3[i]=mesh->cell[st[i]].Ue[k];
                        }
                        wall->URe[k]=weno3R(phi3);
                  }

            }else if(order==5){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi5[i]=mesh->cell[st[i]].Ue[k];
                        }
                        wall->URe[k]=weno5R(phi5);
                  }
            }else if(order==7){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi7[i]=mesh->cell[st[i]].Ue[k];
                        }
                        wall->URe[k]=weno7R(phi7);
                  }
            }else{
                  //order==9
            }



            //LEFT RECONSTRUCTION
            if(wall->nx<TOL4 && wall->nz<TOL4){
                  //y-wall
                  order=wall->cellL->st_sizeY;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellL->stY[j];
                  }
            }else if (wall->nz<TOL4) {
                  //x-wall
                  order=wall->cellL->st_sizeX;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellL->stX[j];
                  }
            }else{
                  //z-wall
                  order=wall->cellL->st_sizeZ;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellL->stZ[j];
                  }
            }

            if(order==1){

                  for(k=0;k<sim->nvar;k++){
                        wall->ULe[k]=wall->cellL->Ue[k];
                  }

            }else if(order==3){

                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi3[i]=mesh->cell[st[i]].Ue[k];
                        }
                        wall->ULe[k]=weno3L(phi3);
                  }

            }else if(order==5){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi5[i]=mesh->cell[st[i]].Ue[k];
                        }
                        wall->ULe[k]=weno5L(phi5);

                  }
            }else if(order==7){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi7[i]=mesh->cell[st[i]].Ue[k];
                        }
                        wall->ULe[k]=weno7L(phi7);
                  }
            }else{
                  //order==7
            }

		uL=wall->ULe[1]/wall->ULe[0];
		uR=wall->URe[1]/wall->URe[0];

		vL=wall->ULe[2]/wall->ULe[0];
		vR=wall->URe[2]/wall->URe[0];

		wL=wall->ULe[3]/wall->ULe[0];
		wR=wall->URe[3]/wall->URe[0];

		wall->pLe=pressure_from_energy(_gamma_, wall->ULe[4], uL, vL, wL, wall->ULe[0], wall->z);
		wall->pRe=pressure_from_energy(_gamma_, wall->URe[4], uR, vR, wR, wall->URe[0], wall->z);

	}

      }

#pragma omp parallel for default(none) private(k,cell) shared(mesh)
	for(i=0;i<mesh->ncells;i++){
		cell=&(mesh->cell[i]);
            if(cell->type!=0){
			cell->S_corr[3] = (cell->w6->pLe-cell->w5->pRe)/cell->dz + _g_*cell->Ue[0];
			//cell->S_corr[3] = (cell->w6->ULe[4]-cell->w5->URe[4])*(_gamma_-1.0)/cell->dz + _g_*cell->Ue[0]; this is only valid for static equilibrium
            }
	}



	return 1;
}


int compute_fluxes(t_mesh *mesh, t_sim *sim){

	double phi3[3],phi5[5],phi7[7]; //auxiliar arrays for the weno reconstruction
	double lambdaMax;
	int order;
	int n,i,j,k;
	int st[9]; //local stencil array
	t_wall *wall;

	mesh->lambda_max=0.0;
	lambdaMax=mesh->lambda_max;
//	wall=mesh->wall;
#pragma omp parallel for default(none) private(wall,phi3,phi5,phi7,order,i,j,k,st) shared(sim,mesh) reduction(max:lambdaMax)
	for(n=0;n<mesh->nwalls;n++){
		wall=&(mesh->wall[n]);

            if(wall->wtype!=0){

            //RIGHT RECONSTRUCTION
            if(wall->nx<TOL4 && wall->nz<TOL4){
                  //y-wall
                  order=wall->cellR->st_sizeY;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellR->stY[j];
                  }
            }else if (wall->nz<TOL4) {
                  //x-wall
                  order=wall->cellR->st_sizeX;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellR->stX[j];
                  }
            }else{
                  //z-wall
                  order=wall->cellR->st_sizeZ;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellR->stZ[j];
                  }
            }
            if(order==1){

                  for(k=0;k<sim->nvar;k++){
                        wall->UR[k]=wall->cellR->U[k];
                  }

            }else if(order==3){

                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi3[i]=mesh->cell[st[i]].U[k];
                        }
                        wall->UR[k]=weno3R(phi3);
                  }

            }else if(order==5){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi5[i]=mesh->cell[st[i]].U[k];
                        }
                        wall->UR[k]=weno5R(phi5);
                  }
            }else if(order==7){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi7[i]=mesh->cell[st[i]].U[k];
                        }
                        wall->UR[k]=weno7R(phi7);
                  }
            }else{
                  //order==9
            }



            //LEFT RECONSTRUCTION
            if(wall->nx<TOL4 && wall->nz<TOL4){
                  //y-wall
                  order=wall->cellL->st_sizeY;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellL->stY[j];
                  }
            }else if (wall->nz<TOL4) {
                  //x-wall
                  order=wall->cellL->st_sizeX;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellL->stX[j];
                  }
            }else{
                  //z-wall
                  order=wall->cellL->st_sizeZ;
                  for(j=0;j<9;j++){
                        st[j]=wall->cellL->stZ[j];
                  }
            }

            if(order==1){

                  for(k=0;k<sim->nvar;k++){
                        wall->UL[k]=wall->cellL->U[k];
                  }

            }else if(order==3){

                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi3[i]=mesh->cell[st[i]].U[k];
                        }
                        wall->UL[k]=weno3L(phi3);
                  }

            }else if(order==5){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi5[i]=mesh->cell[st[i]].U[k];
                        }
                        wall->UL[k]=weno5L(phi5);
                  }
            }else if(order==7){
                  for(k=0;k<sim->nvar;k++){
                        for(i=0;i<order;i++){
                              phi7[i]=mesh->cell[st[i]].U[k];
                        }
                        wall->UL[k]=weno7L(phi7);
                  }
            }else{
                  //order==9
            }


            if(wall->wtype==1){

            //This is to compute fn at each edge
		#if EQUATION_SYSTEM == 2
			#if SOLVER == 0
				compute_euler_HLLE(wall,&lambdaMax);
			#elif SOLVER == 1
				compute_euler_HLLC(wall,&lambdaMax);
			#else 
				compute_euler_HLLS(wall,&lambdaMax,sim);
			#endif
		#elif EQUATION_SYSTEM == 1
                  compute_burgers_flux(wall,&lambdaMax);
            #else
                  compute_linear_flux(wall,&lambdaMax);
            #endif


            }else if(wall->wtype==3){
                  if(wall->boundId==999){
                        printf("%s boundId has been assigned 999, please check. The program will close when pressing a key. \n",ERR);
                        getchar();
                        exit(1);
                  }
            #if EQUATION_SYSTEM == 2
				compute_transmissive_euler(wall,wall->boundId);
            #elif EQUATION_SYSTEM == 1
				compute_burgers_flux(wall,&lambdaMax);
			#else
				compute_linear_flux(wall,&lambdaMax);
			#endif
            }else if(wall->wtype==4){
                  if(wall->boundId==999){
                        printf("%s boundId has been assigned 999, please check. The program will close when pressing a key. \n",ERR);
                        getchar();
                        exit(1);
                  }
            #if EQUATION_SYSTEM == 2
				compute_solid_euler_hlle(wall,&lambdaMax,wall->boundId);
            #endif
            }
            #if EQUATION_SYSTEM == 2
		compute_transport(wall);
            #endif


            }

	}

	mesh->lambda_max=lambdaMax;


	return 1;
}


void compute_transport(t_wall *wall){

	if(wall->fR_star[0]<TOL14){ //negative transport -> information from right hand side
		wall->fR_star[5]=wall->fR_star[0]*wall->UR[5]/wall->UR[0];
		wall->fL_star[5]=wall->fR_star[5];
	}else{ //positive transport -> information from right hand side
		wall->fR_star[5]=wall->fL_star[0]*wall->UL[5]/wall->UL[0];
		wall->fL_star[5]=wall->fR_star[5];
	}


}


void compute_source(t_mesh *mesh){

	int i;
	t_cell *cell;

#pragma omp parallel for default(none) private(cell) shared(mesh)
	for(i=0;i<mesh->ncells;i++){
		cell=&(mesh->cell[i]);
		#if ST==1
            if(cell->type!=0&&cell->st_sizeZ>1){     //This is the implementation of gravity force in -Z direction
			cell->S[3]= -_g_*cell->U[0] + cell->S_corr[3];
			cell->S[4]= -_g_*cell->U[3];
            }
		#elif ST==2
		if(cell->type!=0){     //This is the implementation of gravity force in -Z direction
			cell->S[3]= -_g_*(cell->U[0]-cell->Ue[0]);
			cell->S[4]= -_g_*cell->U[3];
            }
		#else
		if(cell->type!=0){     //This is the implementation of gravity force in -Z direction
			cell->S[3]= -_g_*(cell->U[0]-cell->Ue[0]);
			cell->S[4]= 0.0;
            }
		#endif
	}
}

int update_cell_boundaries(t_mesh *mesh){

	/*
      for(k=0;k<mesh->ncells;k++){

                  p=-1.0/cell[k].dz*(exp(-(cell[k].zc+cell[k].dz/2)) - exp(-(cell[k].zc-cell[k].dz/2)));
			rho=p;

                  if(cell[k].n<1 ){

                        cell[k].U[0]=rho;

                        cell[k].U[4]=p/(_gamma_-1.0)+0.5*rho*(u*u + v*v + w*w); // si la quitamos, error 1.e-12

                  }

                  if(cell[k].n> mesh->zcells-2){


                        cell[k].U[4]=p/(_gamma_-1.0)+0.5*rho*(u*u + v*v + w*w);



      }
	}*/

	return 1;
}


int update_dt(t_mesh *mesh,t_sim *sim){

	double dl;

	dl=MIN(mesh->dx,mesh->dy);
      dl=MIN(dl,mesh->dz);
	sim->dt=sim->CFL*dl/mesh->lambda_max;
	if(sim->dt+sim->t>sim->tf){
		sim->dt=sim->tf-sim->t+TOL14;
	}

	return 1;
}


void mass_calculation(t_mesh *mesh, t_sim *sim){

	int i;
	double massAux;
	double area;

	area=mesh->dx*mesh->dy*mesh->dz;
	massAux=0.0;
#pragma omp parallel for default(none) shared(area,mesh) reduction(+:massAux)
	for(i=0;i<mesh->ncells;i++){
            if(mesh->cell[i].type!=0){
                  massAux+=mesh->cell[i].U[0]*area;
            }
	}
	mesh->mass=massAux;

}

void energy_calculation(t_mesh *mesh, t_sim *sim){

#if EQUATION_SYSTEM==2
	int i;
	double energyAux;
	double area;
	area=mesh->dx*mesh->dy*mesh->dz;
	energyAux=0.0;
#pragma omp parallel for default(none) shared(area,mesh) reduction(+:energyAux)
	for(i=0;i<mesh->ncells;i++){
            if(mesh->cell[i].type!=0){
			#if ST==0||ST==3
			energyAux+=mesh->cell[i].U[4]*area;
			#else
                  energyAux+=(mesh->cell[i].U[4]+mesh->cell[i].U[0]*_g_*mesh->cell[i].zc)*area;
			#endif
            }
	}
	mesh->energy=energyAux;
#else
	mesh->energy=0.0;
#endif
}

void tke_calculation(t_mesh *mesh, t_sim *sim){

	int i;
	double tke_a,u,v,w;
	double volume,volumeT;

	volume=mesh->dx*mesh->dy*mesh->dz;
	volumeT=0.0;
	tke_a=0.0;
	for(i=0;i<mesh->ncells;i++){
            if(mesh->cell[i].type!=0){
                  u=mesh->cell[i].U[1]/mesh->cell[i].U[0];
                  v=mesh->cell[i].U[2]/mesh->cell[i].U[0];
                  w=mesh->cell[i].U[3]/mesh->cell[i].U[0];
                  tke_a+=0.5*mesh->cell[i].U[0]*(u*u + v*v + w*w)*volume;
				  volumeT+=volume;
            }
	}
	mesh->tke=tke_a/volumeT; //average TKE in the domain

}



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
				#if ALLOW_SOLIDS
					update_ghost_cells(sim,mesh,solids);
				#endif
			}else{
				update_cellK1(mesh,sim);
				#if ALLOW_SOLIDS
					update_ghost_cells(sim,mesh,solids);
				#endif                               
			}

		}else if(k==2){
			compute_fluxes(mesh,sim);
			#if ST!=0&&EQUATION_SYSTEM==2
				compute_source(mesh);
			#endif
			update_cellK2(mesh,sim);
			#if ALLOW_SOLIDS
				update_ghost_cells(sim,mesh,solids);
			#endif   				
		}else{ //k=3
			compute_fluxes(mesh,sim);
			#if ST!=0&&EQUATION_SYSTEM==2
				compute_source(mesh);
			#endif
			update_cellK3(mesh,sim);
			#if ALLOW_SOLIDS
				update_ghost_cells(sim,mesh,solids);
			#endif     
		}
	}

}

