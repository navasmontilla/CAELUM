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
  - closures.c

Content:
  -This code file contains all the functions related with the system of equations closures
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
#include "closures.h"



double energy_from_pressure(double gm, double p, double u, double v, double w, double rho, double z){

	double energy;

	#if ST==3
	energy=p/(gm-1.0)+0.5*rho*(u*u + v*v + w*w)+rho*_g_*z;
	#else
	energy=p/(gm-1.0)+0.5*rho*(u*u + v*v + w*w);
	#endif

	return energy;
}

double pressure_from_energy(double gm, double E, double u, double v, double w, double rho, double z){

	double pressure;

	#if ST==3
	pressure=(gm-1.0)*(E-0.5*rho*(u*u + v*v + w*w)-rho*_g_*z);
	#else
	pressure=(gm-1.0)*(E-0.5*rho*(u*u + v*v + w*w));
	#endif


	return pressure;
}
