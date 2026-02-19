/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
- definitions.h

Content:
-This header file contains the macros given to the compiler for the configuration of the software.
 It contains constant values, predefined basic functions and solver options.
 
*/

#ifndef DEFINITIONS_H
  #define DEFINITIONS_H

//Screen messages
#define END "\033[1;32m  =)\033[0m "
#define WAR "\033[1;33m [!]\033[0m "
#define ERR "\033[1;31m [ERROR]\033[0m "
#define OK "\033[1;35m [OK]\033[0m "

//Physical constants
#define PI 3.141592653589793
#define _g_ 9.8
#define _gamma_ 1.4
#define _R_ 287.058
#define _p0_ 1.0E5

//Useful definitions
#define TOL1 1.0E-1
#define TOL2 1.0E-2
#define TOL4 1.0E-4
#define TOL8 1.0E-8
#define TOL14 1.0E-14
#define TOL40 1.0E-40
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define ABS(x) (x < 0 ? -x : x)

//reconstruction method
#define TYPE_REC 0 //This is 0 for WENO, 1 for TENO and 2 for optimal reconstruction
#define _CT_ 1.0e-6
#define epsilon  1.0E-6
#define epsilon2 1.0E-40
#define _Q_ 6.0
#define POSITIVITY 1 //This is 1 if activated, 0 deactivated, positivity fix
#define TOL_RHO 1.0e-8
#define TOL_P 1.0e-8 //2.0e-1 1.0e-8

//Equations
#define EQUATION_SYSTEM 0 // 0: Linear advection, 1: Burgers, 2: Compressible Euler 

//Source terms for Euler
#define ST 0// 0: Source OFF, 1: Source ON (augmented version, needs HLLS), 2: Source ON (perturbation version), 3: Source ON (perturbation version, total energy)

//Multicomponent and multiphase flow
#define MULTICOMPONENT 0 //Activates multicomponent Euler equations (two components with different gamma).
#define MULTI_TYPE 2     //=1 for gamma formulation, =2 for 1/(gamma-1) formulation. ATENTION: Option =2 recommended (see R. Abgrall, S. Karni, Computations of Compressible Multifluids, JCP 169 (2001))

//Solvers
#define SOLVER 1 //0: HLL solver, 1: HLLC solver, 2: HLLS solver

//Debug code
#define DEBUG_MESH 0 //0: no debug; 1: screen info;

//IBM utils
#define ALLOW_SOLIDS 0 //0: no solid cells, 1: STL immersed boundaries, 2: Read list of solid cells, 3: SDF level-set solids
#define _stol_ 3.0 //tolerance for the generation of ghost cell layers. 1.0: 1 layer, 2.0: 2 layers....

//OpenMP configuration
#define NTHREADS 24

//Output files
#define WRITE_VTK 1
#define WRITE_LIST 1
#define WRITE_TKE 0 //write file TKE evolution in time
#define WRITE_PMAX 1

//Printing variables (vtk)
#define print_RHO 1
#define print_VELOCITY 1
#define print_ENERGY 0
#define print_PRESSURE 1
#define print_OVERPRESSURE 0
#define print_SOLUTES 0
#define print_POTENTIALTEM 1

//Reading initial data
#define READ_INITIAL 1



#endif