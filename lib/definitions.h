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
    - definitions.h

  Content:
    -This header file contains the macros given to the compiler
    for the proper configuration of the software.
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
#define TOL4 1.0E-4
#define TOL8 1.0E-8
#define TOL14 1.0E-14
#define TOL40 1.0E-40
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define ABS(x) (x < 0 ? -x : x)

//reconstruction method
#define TYPE_REC 0 //This is 0 for WENO, 1 for TENO and 2 for UWC
#define _CT_ 1.0e-6
#define epsilon  1.0E-6
#define epsilon2 1.0E-40
#define _Q_ 6.0

//Equations
#define EQUATION_SYSTEM 2 // 0: Linear advection, 1: Burgers, 2: Compressible Euler 

//Source terms for Euler
#define ST 3// 0: Source OFF, 1: Source ON (augmented version, needs HLLS), 2: Source ON (perturbation version), 3: Source ON (perturbation version, total energy)

//Multicomponent and multiphase flow
#define MULTICOMPONENT 0 //Activates multicomponent Euler equations (two components with different gamma).
#define MULTI_TYPE 2     //=1 for gamma formulation, =2 for 1/(gamma-1) formulation. ATENTION: Option =2 recommended (see R. Abgrall, S. Karni, Computations of Compressible Multifluids, JCP 169 (2001))

//Solvers
#define SOLVER 0 //0: HLL solver, 1: HLLC solver, 2: HLLS solver

//Debug code
#define DEBUG_MESH 0 //0: no debug; 1: screen info;

//IBM utils
#define ALLOW_SOLIDS 0 //0: no solid cells
#define _stol_ 2.0 //tolerance for the generation of ghost cell layers. 1.0: 1 layer, 2.0: 2 layers....

//OpenMP configuration
#define NTHREADS 32

//Output files
#define WRITE_VTK 1
#define WRITE_LIST 1
#define WRITE_TKE 0 //write file TKE evolution in time

//Printing variables (vtk)
#define print_RHO 0
#define print_VELOCITY 1
#define print_ENERGY 0
#define print_PRESSURE 0
#define print_OVERPRESSURE 1
#define print_SOLUTES 0
#define print_POTENTIALTEM 1

//Reading initial data
#define READ_INITIAL 1



#endif