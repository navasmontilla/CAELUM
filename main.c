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
  - ehow3d.c

Content:
  -This code file contains the main code of the solver.
	It invokes functions which are distributed in libraries within lib folder


*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>

#include "lib/definitions.h"
#include "lib/structures.h"
#include "lib/closures.h"
#include "lib/ibmutils.h"
#include "lib/mathutils.h"
#include "lib/numcore.h"
#include "lib/postproc.h"
#include "lib/preproc.h"
#include "lib/reconst.h"
#include "lib/solvers.h"



int main(int argc, char * argv[]){

	t_mesh *mesh;
	t_sim *sim;
	t_solid *solids;
	char vtkfile[1024];
	char listfile[1024];
	double tf,timeac;
	int nIt;
#if WRITE_TKE
	double timeac2,tTke;
      FILE *file_tke;
#endif


    if (argc < 2) {
        printf("%s A folder path must be passed as follows: %s <folder_path>\n", ERR, argv[0]);
        return 1; // Exit 
    }
    const char *folder_path = argv[1];

#ifdef _OPENMP
	omp_set_num_threads(NTHREADS);
     	printf("The number of threads is set to %d.\n",NTHREADS);
#pragma omp parallel
     printf("Thread %d of %d is checked.\n", omp_get_thread_num(), omp_get_num_threads());

#endif


      ////////////////////////////////////////////////////
	///////////// M E M O R Y  A L L O C. //////////////
	////////////////////////////////////////////////////

	//Mesh allocation
	mesh=(t_mesh*)malloc(sizeof(t_mesh));

	//Simulation allocation
	sim=(t_sim*)malloc(sizeof(t_sim));

	//Solids allocation
	solids=(t_solid*)malloc(sizeof(t_solid));


	////////////////////////////////////////////////////
	////////////// P R E - P R O C E S S ///////////////
	////////////////////////////////////////////////////
	
	read_config(mesh,sim,folder_path);
	print_info(mesh,sim,folder_path);

#if EQUATION_SYSTEM == 2
	sim->nvar=6;
#else
	sim->nvar=1;
#endif

	if(sim->order==1){
		sim->rk_steps=1;
	}else{
		sim->rk_steps=3;
	}

	mesh->dx= mesh->Lx/mesh->xcells;
	mesh->dy= mesh->Ly/mesh->ycells;
      mesh->dz= mesh->Lz/mesh->zcells;

	timeac=0.0;
#if WRITE_TKE
	timeac2=0.0;
      tTke=0.05;
#endif
      

	////////////////////////////////////////////////////
	////////// I N I T I A L I Z A T I O N /////////////
	////////////////////////////////////////////////////
	
	create_mesh(mesh,sim);
#if ALLOW_SOLIDS
	read_solids(mesh,solids,folder_path);
#else
	solids->nsolid=0;
#endif
      assign_cell_type(mesh,solids);
      update_stencils(mesh,sim);
      assign_wall_type(mesh);    
	update_initial(mesh,sim,folder_path);

#if ALLOW_SOLIDS
      assign_image_cells(mesh,solids);
      update_ghost_cells(sim,mesh,solids);
      update_wall_type(mesh,solids);
      printf("%s Image points have been defined and ghost cell values have been computed \n",OK);
#endif

#if EQUATION_SYSTEM == 0
	set_velocity(mesh,sim);
#endif
//	//update_boundaries(mesh,bc); //??
	
	snprintf(vtkfile, sizeof(vtkfile),"%s/out/inital_geo_mesh.vtk", folder_path);
	write_geo_vtk(mesh,vtkfile);
	snprintf(vtkfile, sizeof(vtkfile),"%s/out/state000.vtk", folder_path);
	write_vtk(mesh,vtkfile);
	snprintf(listfile, sizeof(listfile),"%s/out/state000.out", folder_path);
	write_list(mesh,listfile);
	snprintf(listfile, sizeof(listfile),"%s/out/list_eq.out", folder_path);
	write_list_eq(mesh,listfile);
      printf("\n");
      printf(" T= 0.0e+0. Initial data printed. Starting time loop.\n");


#if WRITE_TKE == 1
	snprintf(listfile, sizeof(listfile),"%s/out/tke.out", folder_path);
	file_tke=fopen(listfile,"w");
#endif


	mass_calculation(mesh,sim);
	energy_calculation(mesh,sim);
	mesh->mass0=mesh->mass;
	mesh->energy0=mesh->energy;

	////////////////////////////////////////////////////
	////////////// C A L C U L A T I O N ///////////////
	////////////////////////////////////////////////////

	tf=sim->tf;
      sim->t=0.0;
	nIt=0;

	#if ST!=0&&EQUATION_SYSTEM==2
      equilibrium_reconstruction(mesh,sim);
	#endif

	while(sim->t<tf){
		
		update_solution(mesh,sim,solids,sim->rk_steps); //updates all variables one dt

		if(mesh->cell_bc_flag!=1){
			update_cell_boundaries(mesh);
		}

		////////////////////////////////////////////////////
		////////////////// P R I N T I N G /////////////////
		////////////////////////////////////////////////////

            timeac=timeac+sim->dt;
		sim->t+=sim->dt;	//Time for the next time step

		if(timeac>sim->tVolc){
			
			screen_info(mesh,sim);
			
			#if WRITE_VTK
			snprintf(vtkfile, sizeof(vtkfile), "%s/out/state%03d.vtk", folder_path, nIt + 1);
			write_vtk(mesh,vtkfile);
			#endif
			#if WRITE_LIST
			snprintf(listfile, sizeof(listfile), "%s/out/state%03d.out", folder_path, nIt + 1);
			write_list(mesh,listfile);
			#endif
			
			nIt++;		//Number of iteration for the next time step
			timeac=0.0;
		}

		#if WRITE_TKE
            timeac2=timeac2+sim->dt;
		if(timeac2>tTke){
            tke_calculation(mesh,sim);
            fprintf(file_tke,"%14.14e %14.14e\n",sim->t,mesh->tke);
			timeac2=0.0;
		}
		#endif


	}

      printf(" \n");
      printf(" Final time is T= %14.14e \n \n",sim->t);

      if(timeac>TOL14){
	snprintf(vtkfile, sizeof(vtkfile), "%s/out/state%03d.vtk", folder_path, nIt + 1);
      write_vtk(mesh,vtkfile);
	snprintf(listfile, sizeof(listfile), "%s/out/state%03d.out", folder_path, nIt + 1);
	write_list(mesh,listfile);
      }

      #if WRITE_TKE
      fclose(file_tke);
      #endif



	printf("\n%s Simulation completed!\n",END);



	return 1;

}
