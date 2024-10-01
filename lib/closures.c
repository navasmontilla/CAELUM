/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - closures.c

Content:
  -This code file contains all the functions related with the closure equations for the pressure
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
